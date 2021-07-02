#########################################################
# CancerCLAS Multiclass Prediction
# Empirical Data
# Naive No Variation, Naive Bootstrap, Weighted Bootstrap
# Results Compiling / Writing Output
#########################################################
options(scipen = 999)

library(MASS)
library(tidyverse)
library(glmnet)
library(nnet)
library(gbm)
library(mgcv) 
library(xgboost)
library(polspline)
library(randomForest)
library(caret)
library(pROC)
library(survival)
library(doParallel)
library(ranger)
detectCores()
library(nnls)
library(SuperLearner)
library(survminer)
library(ggfortify)
library(rms)
library(pec)
library(riskRegression)
library(xtable)


# source misc_funs.R: contains function for calculating class measures
source("~/misc_funs_empirical_v3.R")

# seed and read in
set.seed(33)
readpath <- c("~/read_path/")

data_dev_thlds <- readRDS(paste0(readpath, "thlds",".RDS"))
val_label_results <- readRDS(paste0(readpath, "val_wtd_cov_amb_in_sample", ".RDS"))
val_naive_sample_coverage <- readRDS(paste0(readpath, "naive_sample_coverage", ".RDS"))
val_naive_pred_in_sample <- readRDS(paste0(readpath, "naive_class_perf_in_sample", ".RDS"))
sample_km_obs <- readRDS(paste0(readpath, "km_obs_in_sample",".RDS"))
sample_km_pred_naive <- readRDS(paste0(readpath, "km_naive_pred_in_sample",".RDS"))
bs_measures_list <- readRDS(paste0(readpath, "bs_measures_list", ".RDS"))

# write tables to .txt files
wrFun <- function(df, name){
  write.table(df, file=paste0(readpath, name, ".txt"),sep=",", quote=F, row.names=F)
}

n_val <- length(val_label_results$Hstar$mn$g1)
# true class counts
n_g1 <- 2169; n_g2 <- 4851; n_g3 <- 7599
alg_names <- c("mn","lasso","ridge","enet","rf","gam","xgb")
# thlds
tmp_thlds <- round(do.call(rbind.data.frame, data_dev_thlds),3)
colnames(tmp_thlds) <- c("g1","g2","g3")
tmp_thlds <- tmp_thlds %>% mutate(alg=alg_names) %>% select(alg, everything()) %>%
  arrange(alg)
# Table C14: Weighted labeling thresholds.
print(xtable(tmp_thlds), include.rownames=F)
wrFun(tmp_thlds, "/thlds"); rm(tmp_thlds)

# wtd amb total in sample 
# counts
tmp_wtd_amb <- val_label_results$amb_total
tmp_wtd_amb <- do.call(rbind.data.frame, tmp_wtd_amb)
tmp_wtd_amb <- tmp_wtd_amb %>% mutate(alg=alg_names) %>% select(alg, everything())
# proportions
tmp_wtd_amb_pct <- tmp_wtd_amb %>% mutate(l0=round(l0/n_val,2),
                                          l1=round(l1/n_val,2),
                                          l2=round(l2/n_val,2),
                                          l3=round(l3/n_val,2))
wrFun(tmp_wtd_amb, "/val_wtd_amb_total_in_sample"); wrFun(tmp_wtd_amb_pct, "/val_wtd_amb_total_pct_in_sample"); rm(tmp_wtd_amb)
# wtd amb by true class in sample
tmp_wtd_amb <- val_label_results$amb_class
tmp_wtd_amb_tot <- do.call(rbind.data.frame, tmp_wtd_amb)
tmp_wtd_amb_tot <- cbind(alg_class=rownames(tmp_wtd_amb_tot),data.frame(tmp_wtd_amb_tot,row.names=NULL))
# proportions
acFun <- function(data){
  tmp <- data %>% mutate(class=c(1,2,3))
  tmp_g1 <- tmp %>% filter(class==1) %>% mutate_at(vars(-class), funs(./n_g1))
  tmp_g2 <- tmp %>% filter(class==2) %>% mutate_at(vars(-class), funs(./n_g2))
  tmp_g3 <- tmp %>% filter(class==3) %>% mutate_at(vars(-class), funs(./n_g3))
  tmp <- rbind.data.frame(tmp_g1, tmp_g2, tmp_g3)
  tmp <- tmp %>% select(class, everything())
  tmp <- round(tmp,2)
  return(tmp)
}
tmp_wtd_amb_pct <- lapply(tmp_wtd_amb, acFun)
tmp_wtd_amb_pct <- do.call(rbind.data.frame, tmp_wtd_amb_pct)
tmp_wtd_amb_pct <- cbind(alg_class=rownames(tmp_wtd_amb_pct),data.frame(tmp_wtd_amb_pct,row.names=NULL))
wrFun(tmp_wtd_amb_tot, "/val_wtd_amb_class_in_sample"); wrFun(tmp_wtd_amb_pct, "/val_wtd_amb_class_pct_in_sample"); rm(tmp_wtd_amb)
# naive coverage in sample -- use bs version instead
tmp_naive_cov <- round(do.call(rbind.data.frame, val_naive_sample_coverage),3)
colnames(tmp_naive_cov) <- c("g1","g2","g3")
tmp_naive_cov <- tmp_naive_cov %>% mutate(alg=alg_names) %>% select(alg, everything())
wrFun(tmp_naive_cov, "/val_naive_cov_in_sample"); rm(tmp_naive_cov)
# naive calibration
calFun <- function(measure){
  tmp <- lapply(val_naive_pred_in_sample,`[[`,measure)
  tmp <- lapply(tmp, function(x){do.call(cbind.data.frame,x)})
  tmp <- do.call(rbind.data.frame, tmp)
  tmp <- tmp  %>% mutate(alg=c(rep(alg_names[1],30), rep(alg_names[2],30),
                               rep(alg_names[3],30), rep(alg_names[4],30),
                               rep(alg_names[5],30), rep(alg_names[6],30),
                               rep(alg_names[7],30))) %>% select(alg, everything())
  return(tmp)
}
tmp_cal <- calFun("cal")
wrFun(tmp_cal, "/val_naive_in_sample_cal");rm(tmp_cal)
# naive class performance in sample
cmFun <- function(measure){
  tmp <- lapply(val_naive_pred_in_sample,`[[`,measure)
  tmp <- lapply(tmp, function(x){do.call(cbind.data.frame,x)})
  tmp <- round(do.call(rbind.data.frame, tmp),2)
  tmp <- tmp %>% mutate(alg=alg_names) %>% select(alg, everything())
  return(tmp)
}
tmp_acc <- cmFun("acc")
# Table C15: Classification accuracy.
print(xtable(arrange(tmp_acc, alg)),include.rownames=F)
wrFun(tmp_acc, "/val_naive_in_sample_acc"); rm(tmp_acc)
tmp_sens <- cmFun("sens")
# Table C16: Classification sensitivity
print(xtable(arrange(tmp_sens, alg)),include.rownames=F)
wrFun(tmp_sens, "/val_naive_in_sample_sens"); rm(tmp_sens)
tmp_spec <- cmFun("spec")
# Table C17: Classification specificity
print(xtable(arrange(tmp_spec, alg)),include.rownames=F)
wrFun(tmp_spec, "/val_naive_in_sample_spec"); rm(tmp_spec)
tmp_ppv <- cmFun("ppv")
# Table C18: Classification PPV
print(xtable(arrange(tmp_ppv, alg)),include.rownames=F)
wrFun(tmp_ppv, "/val_naive_in_sample_ppv"); rm(tmp_ppv)

# KM based on observed class in sample
# median
tmp_km_med <- sample_km_obs$sample_median_obs
tmp_km_med <- cbind.data.frame(tmp_km_med$quantile, tmp_km_med$lower, tmp_km_med$upper)
colnames(tmp_km_med) <- c("median","lb","ub")
tmp_km_med <- tmp_km_med %>% mutate(class=c("c1","c2","c3")) %>% 
  select(class, everything()) %>% mutate(ci_width=ub-lb)
wrFun(tmp_km_med, "/km_obs_in_sample_median"); rm(tmp_km_med)
# times
tmp_km_times <- sample_km_obs$sample_times_obs
# d90
tmp_km_d90 <- round(cbind.data.frame(tmp_km_times$d90$surv, tmp_km_times$d90$lb, tmp_km_times$d90$ub),2)
colnames(tmp_km_d90) <- c("mean","lb","ub")
tmp_km_d90 <- tmp_km_d90 %>% mutate(class=c("c1","c2","c3"), time=rep("90",3)) %>% 
  select(time,class, everything()) %>% mutate(ci_width=ub-lb)
# d365
tmp_km_d365 <- round(cbind.data.frame(tmp_km_times$d365$surv, tmp_km_times$d365$lb, tmp_km_times$d365$ub),2)
colnames(tmp_km_d365) <- c("mean","lb","ub")
tmp_km_d365 <- tmp_km_d365 %>% mutate(class=c("c1","c2","c3"),time=rep("365",3)) %>% 
  select(time,class, everything()) %>% mutate(ci_width=ub-lb)

wrFun(tmp_km_d90, "/km_obs_in_sample_d90"); wrFun(tmp_km_d365, "/km_obs_in_sample_d365")
rm(tmp_km_d90, tmp_km_d365, tmp_km_times)

# KM based on naive predicted class in sample
kmMedFun <- function(data){
  tmp_km_med <- data[["sample_median_pred"]]
  tmp_median <- as.data.frame.matrix(t(tmp_km_med$quantile))
  tmp_lb <- as.data.frame.matrix(t(tmp_km_med$lower))
  tmp_ub <- as.data.frame.matrix(t(tmp_km_med$upper))
  colnames(tmp_median) <- c("g1","g2","g3")
  colnames(tmp_lb) <- c("lb_g1","lb_g2","lb_g3")
  colnames(tmp_ub) <- c("ub_g1","ub_g2","ub_g3")
  tmp_km_med <- cbind.data.frame(tmp_median, tmp_lb, tmp_ub)
  tmp_km_med <- tmp_km_med %>% mutate(w_g1=ub_g1-lb_g1,w_g2=ub_g2-lb_g2,w_g3=ub_g3-lb_g3)
  tmp_km_med <- tmp_km_med %>% select(g1,lb_g1,ub_g1,w_g1,g2,lb_g2,ub_g2,w_g2,g3,lb_g3,ub_g3,w_g3)
  return(tmp_km_med)
}
tmp_km_med <- lapply(sample_km_pred_naive, kmMedFun)
tmp_km_med <- do.call(rbind.data.frame, tmp_km_med)
tmp_km_med <- cbind(alg=rownames(tmp_km_med),data.frame(tmp_km_med,row.names=NULL))
wrFun(tmp_km_med, "/km_naive_pred_in_sample_median")

# times
kmTimesFun <- function(data,time){
  tmp <- data[["sample_times_pred"]]
  tmp <- tmp[[time]]
  tmp <- cbind.data.frame(t(tmp$surv), t(tmp$lb), t(tmp$ub))
  tmp <- round(tmp,2)
  colnames(tmp) <- c("g1","g2","g3","lb_g1","lb_g2","lb_g3","ub_g1","ub_g2","ub_g3")
  tmp <- tmp %>% mutate(w_g1=ub_g1-lb_g1,w_g2=ub_g2-lb_g2,w_g3=ub_g3-lb_g3)
  tmp <- tmp %>% select(g1,lb_g1,ub_g1,w_g1,g2,lb_g2,ub_g2,w_g2,g3,lb_g3,ub_g3,w_g3)
  return(tmp)
}
# d90
tmp_km_d90 <- lapply(sample_km_pred_naive, kmTimesFun, "d90")
tmp_km_d90 <- do.call(rbind.data.frame, tmp_km_d90)
tmp_km_d90 <- cbind(alg=rownames(tmp_km_d90),data.frame(tmp_km_d90,row.names=NULL))
wrFun(tmp_km_d90, "/km_naive_pred_in_sample_d90")
# d365
tmp_km_d365 <- lapply(sample_km_pred_naive, kmTimesFun, "d365")
tmp_km_d365 <- do.call(rbind.data.frame, tmp_km_d365)
tmp_km_d365 <- cbind(alg=rownames(tmp_km_d365),data.frame(tmp_km_d365,row.names=NULL))
wrFun(tmp_km_d365, "/km_naive_pred_in_sample_d365")

### BS-based measures ###
# naive class coverage
tmp_cc <- lapply(bs_measures_list, `[[`, "class_coverage_naive")
tmp_cc <- do.call(rbind.data.frame, tmp_cc)
tmp_cc <- apply(tmp_cc, 2, round, digits=2)
tmp_cc<- cbind(alg=rownames(tmp_cc),data.frame(tmp_cc,row.names=NULL))
# for latex think these will look beter with bounds as add'l rows
wrFun(tmp_cc, "/bs_naive_class_coverage")

boundsFun <- function(data){
  tmp <- data %>% select(alg,starts_with("lb"), starts_with("ub"))
  tmp <- tmp %>% mutate(mean_g1=paste0("(",lb_g1, ", ",ub_g1,")"),
                                mean_g2=paste0("(",lb_g2, ", ",ub_g2,")"),
                                mean_g3=paste0("(",lb_g3, ", ",ub_g3,")")) %>%
    select(alg, mean_g1, mean_g2, mean_g3)
  out <- data %>% select(alg, mean_g1, mean_g2, mean_g3) %>% rbind(., tmp) %>%
    arrange(alg)
  return(out)
}
tmp <- boundsFun(tmp_cc)
# Table C13: Bootstrap-based coverage (naive portion)
print(xtable(tmp),include.rownames=F); rm(tmp, tmp_cc)
# need to check if these means are across BS or taken from sample
# wtd class coverage
tmp_cc <- do.call(rbind.data.frame, tmp_cc)
tmp_cc <- apply(tmp_cc,2,round, digits=2)
tmp_cc<- cbind(alg=rownames(tmp_cc),data.frame(tmp_cc,row.names=NULL))
wrFun(tmp_cc, "/bs_wtd_class_coverage")
tmp <- boundsFun(tmp_cc)
# Table C13: Bootstrap-based coverage (Wtd portion)
print(xtable(tmp), include.rownames=F); rm(tmp, tmp_cc)

# ambiguity
tmp_amb <- lapply(bs_measures_list, `[[`, "ambiguity")
aFun <- function(meas){
  tmp <- lapply(tmp_amb, `[[`, meas)
  tmp <- round(do.call(rbind.data.frame, tmp))
  tmp <- cbind(alg=rownames(tmp),data.frame(tmp,row.names=NULL))
  return(tmp)
}
# total
tmp_amb_all <- aFun("amb_all")
# pct
tmp_amb_all_pct <- tmp_amb_all %>% mutate_at(vars(-alg), funs(./n_val)) %>% mutate_at(vars(-alg), ~round(.,2))
wrFun(tmp_amb_all, "/bs_wtd_amb_total"); wrFun(tmp_amb_all_pct, "/bs_wtd_amb_total_pct")
boundsFun2 <- function(data){
  tmp <- data %>% select(alg,starts_with("lb"), starts_with("ub"))
  tmp <- tmp %>% mutate(mean_l0=paste0("(",lb_l0, ", ",ub_l0,")"),
                        mean_l1=paste0("(",lb_l1, ", ",ub_l1,")"),
                        mean_l2=paste0("(",lb_l2, ", ",ub_l2,")"),
                        mean_l3=paste0("(",lb_l3, ", ",ub_l3,")"))%>%
    select(alg, mean_l0, mean_l1, mean_l2, mean_l3)
  out <- data %>% select(alg,mean_l0, mean_l1, mean_l2, mean_l3) %>% rbind(., tmp) %>%
    arrange(alg)
  return(out)
}
tmp <- boundsFun2(tmp_amb_all_pct)
rm(tmp_amb_all,tmp_amb_all_pct, tmp)
# by true class 
# class 1
tmp_g1 <- aFun("amb_g1")
# pct -- for pct just do the means (using sum of means for each class, which is same as emp sample)
tmp_g1_pct <- tmp_g1 %>% 
  mutate_at(vars(-alg),funs(./n_g1)) %>% mutate_at(vars(-alg), funs(round(.,2)))
wrFun(tmp_g1, "/bs_wtd_amb_g1"); wrFun(tmp_g1_pct, "/bs_wtd_amb_g1_pct"); 
tmp_g1_pct <- tmp_g1_pct %>% select(alg, starts_with("mean"))
rm(tmp_g1, tmp_g1_pct)
# class 2
tmp_g2 <- aFun("amb_g2")
tmp_g2_pct <- tmp_g2 %>% 
  mutate_at(vars(-alg), funs(./n_g2)) %>% mutate_at(vars(-alg),funs(round(.,2)))
wrFun(tmp_g2, "/bs_wtd_amb_g2"); wrFun(tmp_g2_pct, "/bs_wtd_amb_g2_pct")
tmp_g2_pct <- tmp_g2_pct %>% select(alg, starts_with("mean"))
rm(tmp_g2, tmp_g2_pct)
# class 3
tmp_g3 <- aFun("amb_g3")
tmp_g3_pct <- tmp_g3 %>% 
  mutate_at(vars(-alg), funs(./n_g3)) %>% mutate_at(vars(-alg),funs(round(.,2)))
wrFun(tmp_g3, "/bs_wtd_amb_g3"); wrFun(tmp_g3_pct, "/bs_wtd_amb_g3_pct")
tmp_g3_pct <- tmp_g3_pct %>% select(alg, starts_with("mean"))
rm(tmp_g3, tmp_g3_pct)

# BS calibration
calBSFun <- function(method,measure){
  tmp <- lapply(bs_measures_list,`[[`,method)
  tmp <- lapply(tmp, `[[`,measure)
  tmp <- lapply(tmp, function(x){do.call(cbind.data.frame,x)})
  tmp <- do.call(rbind.data.frame, tmp)
  tmp <- tmp  %>% mutate(alg=c(rep(alg_names[1],30), rep(alg_names[2],30),
                               rep(alg_names[3],30), rep(alg_names[4],30),
                               rep(alg_names[5],30), rep(alg_names[6],30),
                               rep(alg_names[7],30))) %>% select(alg, everything())
  return(tmp)
}
tmp_cal <- calBSFun("class_measures_naive_pct","cal")
wrFun(tmp_cal, "/bs_naive_pct_cal")
tmp_cal_wtd <- calBSFun("class_measures_wtd","cal")
wrFun(tmp_cal_wtd, "/bs_wtd_cal")

# naive (pct) class measures
# naive class performance in sample
cmBSFun <- function(method,measure){
  tmp <- lapply(bs_measures_list,`[[`,method)
  tmp <- lapply(tmp, `[[`,measure)
  tmp <- lapply(tmp, function(x){do.call(cbind.data.frame,x)})
  tmp <- round(do.call(rbind.data.frame, tmp),2)
  tmp <- tmp %>% mutate(alg=alg_names) %>% select(alg, everything())
  return(tmp)
}

boundsFun3 <- function(data){
    tmp <- data %>% select(alg,starts_with("lb"), starts_with("ub"))
    tmp <- tmp %>% mutate(all=paste0("(",lb_all, ", ",ub_all,")"),
                          g1=paste0("(",lb_g1, ", ",ub_g1,")"),
                          g2=paste0("(",lb_g2, ", ",ub_g2,")"),
                          g3=paste0("(",lb_g3, ", ",ub_g3,")"))%>%
      select(alg, all, g1, g2, g3)
    out <- data %>% select(alg,all, g1, g2, g3) %>% rbind(., tmp) %>%
      arrange(alg)
    return(out)
  }
tmp_acc <- cmBSFun("class_measures_naive_pct","acc")
wrFun(tmp_acc, "/bs_naive_pct_acc")
tmp_acc <- boundsFun3(tmp_acc)
# Table C15: Classification Accuracy
print(xtable(arrange(tmp_acc, alg)),include.rownames=F); rm(tmp_acc)

tmp_sens <- cmBSFun("class_measures_naive_pct","sens")
wrFun(tmp_sens, "/bs_naive_pct_sens")
tmp_sens <- boundsFun3(tmp_sens)
# Table C16: Classification sensitivity
print(xtable(arrange(tmp_sens, alg)),include.rownames=F); rm(tmp_sens)

tmp_spec <- cmBSFun("class_measures_naive_pct","spec")
wrFun(tmp_spec, "/bs_naive_pct_spec")
tmp_spec <- boundsFun3(tmp_spec)
# Table C17: Classification specificity
print(xtable(arrange(tmp_spec, alg)),include.rownames=F); rm(tmp_spec)

tmp_ppv <- cmBSFun("class_measures_naive_pct","ppv")
wrFun(tmp_ppv, "/bs_naive_pct_ppv")
tmp_ppv <- boundsFun3(tmp_ppv)
# Table C18: Classification ppv
print(xtable(arrange(tmp_ppv, alg)),include.rownames=F); rm(tmp_ppv)

ccBSFun <- function(method,grp){
  tmp <- lapply(bs_measures_list,`[[`,method)
  tmp <- lapply(tmp,`[[`,grp)
  tmp <- round(do.call(rbind.data.frame, tmp))
  tmp <- tmp %>% mutate(alg=alg_names) %>% select(alg, everything())
  return(tmp)
}
tmp_tp <- ccBSFun("class_measures_naive_pct","tp")
wrFun(tmp_tp, "/bs_naive_pct_tp"); rm(tmp_tp)
tmp_tn <- ccBSFun("class_measures_naive_pct","tn")
wrFun(tmp_tn, "/bs_naive_pct_tn"); rm(tmp_tn)
tmp_fp <- ccBSFun("class_measures_naive_pct","fp")
wrFun(tmp_fp, "/bs_naive_pct_fp"); rm(tmp_fp)
tmp_fn <- ccBSFun("class_measures_naive_pct","fn")
wrFun(tmp_fn, "/bs_naive_pct_fn"); rm(tmp_fn)
# wtd class measures
tmp_acc <- cmBSFun("class_measures_wtd","acc")
wrFun(tmp_acc, "/bs_wtd_acc")
tmp_acc <- boundsFun3(tmp_acc)
# Table C15: classification accuracy
print(xtable(arrange(tmp_acc, alg)),include.rownames=F); rm(tmp_acc)

tmp_sens <- cmBSFun("class_measures_wtd","sens")
wrFun(tmp_sens, "/bs_wtd_sens")
tmp_sens <- boundsFun3(tmp_sens)
# Table C16: classification sensitivity
print(xtable(arrange(tmp_sens, alg)),include.rownames=F); rm(tmp_sens)

tmp_spec <- cmBSFun("class_measures_wtd","spec")
wrFun(tmp_spec, "/bs_wtd_spec")
tmp_spec <- boundsFun3(tmp_spec)
# Table C17: Classification specificity
print(xtable(arrange(tmp_spec, alg)),include.rownames=F); rm(tmp_spec)

tmp_ppv <- cmBSFun("class_measures_wtd","ppv")
wrFun(tmp_ppv, "/bs_wtd_ppv")
tmp_ppv <- boundsFun3(tmp_ppv)
# Table C18: Classficiation ppv
print(xtable(arrange(tmp_ppv, alg)),include.rownames=F); rm(tmp_ppv)


tmp_tp <- ccBSFun("class_measures_wtd","tp")
wrFun(tmp_tp, "/bs_wtd_tp"); rm(tmp_tp)
tmp_tn <- ccBSFun("class_measures_wtd","tn")
wrFun(tmp_tn, "/bs_wtd_tn"); rm(tmp_tn)
tmp_fp <- ccBSFun("class_measures_wtd","fp")
wrFun(tmp_fp, "/bs_wtd_fp"); rm(tmp_fp)
tmp_fn <- ccBSFun("class_measures_wtd","fn")
wrFun(tmp_fn, "/bs_wtd_fn"); rm(tmp_fn)
# km naive (pct)

# KM based on naive predicted class 
kmBSFun <- function(method, measure){
  tmp <- lapply(bs_measures_list, `[[`,method)
  tmp <- lapply(tmp,`[[`,measure)
  tmp <- do.call(rbind.data.frame, tmp)
  tmp <- cbind(alg=rownames(tmp),data.frame(tmp,row.names=NULL))
  tmp <- tmp %>% mutate(w_g1=ub_g1-lb_g1,w_g2=ub_g2-lb_g2,w_g3=ub_g3-lb_g3)
  tmp <- tmp %>% select(alg, g1,lb_g1,ub_g1,w_g1,g2,lb_g2,ub_g2,w_g2,g3,lb_g3,ub_g3,w_g3)
  tmp <- tmp %>% mutate_at(vars(-alg), funs(round(.,2)))
  return(tmp)
}
tmp_naive_med <- kmBSFun(method="km_naive_pct",measure="median")
tmp_naive_med <- tmp_naive_med %>% mutate_at(vars(-alg),funs(round(.)))
tmp_naive_d90 <- kmBSFun(method="km_naive_pct",measure="d90")
tmp_naive_d365 <- kmBSFun(method="km_naive_pct",measure="d365")

wrFun(tmp_naive_med, "/bs_naive_pct_km_median")
wrFun(tmp_naive_d90, "/bs_naive_pct_km_d90")
wrFun(tmp_naive_d365, "/bs_naive_pct_km_d365")


tmp_naive_med <- boundsFun4(tmp_naive_med)
print(xtable(arrange(tmp_naive_med, alg)), include.rownames=F); rm(tmp_naive_med)
tmp_naive_d90 <- boundsFun4(tmp_naive_d90)
print(xtable(arrange(tmp_naive_d90, alg)), include.rownames=F); rm(tmp_naive_d90)
tmp_naive_d365 <- boundsFun4(tmp_naive_d365)
print(xtable(arrange(tmp_naive_d365, alg)), include.rownames=F); rm(tmp_naive_d365)

# km wtd
tmp_wtd_med <- kmBSFun(method="km_wtd",measure="median")
tmp_wtd_med <- tmp_wtd_med %>% mutate_at(vars(-alg),funs(round(.)))
tmp_wtd_d90 <- kmBSFun(method="km_wtd",measure="d90")
tmp_wtd_d365 <- kmBSFun(method="km_wtd",measure="d365")
wrFun(tmp_wtd_med, "/bs_wtd_km_median")
wrFun(tmp_wtd_d90, "/bs_wtd_km_d90")
wrFun(tmp_wtd_d365, "/bs_wtd_km_d365")


tmp_wtd_med <- boundsFun4(tmp_wtd_med)
print(xtable(arrange(tmp_wtd_med, alg)), include.rownames=F); rm(tmp_wtd_med)
tmp_wtd_d90 <- boundsFun4(tmp_wtd_d90)
print(xtable(arrange(tmp_wtd_d90, alg)), include.rownames=F); rm(tmp_wtd_d90)
tmp_wtd_d365 <- boundsFun4(tmp_wtd_d365)
print(xtable(arrange(tmp_wtd_d365, alg)), include.rownames=F); rm(tmp_wtd_d365)


# bias -- comparing naive no-boot, naive boot, and wtd boot to observed
# for bs: data is bs, method is km_wtd or naive_pct, measures are d90 etc
# for in sample, data is sample_km_pred_naive, method is sample_times or sample_median
# and measure is d90 etc
biasFun <- function(data,method, measure, obs_measure){
  tmp <- lapply(data, `[[`,method)
  tmp <- lapply(tmp,`[[`,measure)
  if(method %in% c("km_wtd","km_naive_pct")){
    tmp <- lapply(tmp,`[`,c("g1","g2","g3"))
  }else if(measure %in% c("quantile","lower","upper")) {
    tmp <- lapply(tmp, function(x){as.data.frame.matrix(t(x))})
  }else{
    tmp <- lapply(tmp, `[[`,"surv")
  }
  tmp_bias <- lapply(tmp, function(x){x-obs_measure})
  tmp_bias <- round(do.call(rbind.data.frame, tmp_bias),3)
  colnames(tmp_bias) <- c("g1","g2","g3")
  tmp_bias <- tmp_bias %>% mutate(alg=alg_names) %>% select(alg, everything())
  return(tmp_bias)
}
pctBiasFun <- function(data, method, measure, obs_measure){
  tmp <- lapply(data, `[[`, method)
  tmp <- lapply(tmp, `[[`, measure)
  if(method %in% c("km_wtd","km_naive_pct")){
    tmp <- lapply(tmp,`[`,c("g1","g2","g3"))
  }else if(measure %in% c("quantile","lower","upper")) {
    tmp <- lapply(tmp, function(x){as.data.frame.matrix(t(x))})
  }else{
    tmp <- lapply(tmp, `[[`,"surv")
  }
  pct_bias <- lapply(tmp, function(x){100*((x-obs_measure)/obs_measure)})
  pct_bias <- round(do.call(rbind.data.frame, pct_bias),1)
  colnames(pct_bias) <- c("g1","g2","g3")
  pct_bias <- pct_bias %>% mutate(alg=alg_names) %>% select(alg, everything())
  return(pct_bias)
}
# d90
# calcualte naive bias here -- bootstrap versions were done in bootstrap
# clean bs versions
bsBiasFun <- function(method, measure, grp){
  tmp <- lapply(bs_measures_list,`[[`,method)
  tmp <- lapply(tmp, `[[`, measure)
  tmp <- do.call(rbind.data.frame, tmp)
  tmp <- tmp %>% mutate(alg=alg_names, setting=grp)
  return(tmp)
}
d90_obs <- sample_km_obs$sample_times_obs$d90$surv
naive_sample_d90_bias <- biasFun(sample_km_pred_naive,"sample_times_pred","d90",d90_obs)
wrFun(naive_sample_d90_bias, "/naive_sample_d90_bias")
print(xtable(arrange(naive_sample_d90_bias,alg), digits=3), include.rownames=F); rm(naive_sample_d90_bias)

bs_naive_pct_d90_bias <- bsBiasFun("km_naive_pct","bias_90","naive_bs_pct")
wrFun(bs_naive_pct_d90_bias, "/bs_naive_pct_d90_bias")
bs_wtd_d90_bias <- bsBiasFun("km_wtd","bias_90","bs_wtd")
wrFun(bs_wtd_d90_bias, "/bs_wtd_d90_bias")
# d90 pct bias
naive_sample_d90_bias_pct <- pctBiasFun(sample_km_pred_naive,"sample_times_pred","d90",d90_obs)
bs_naive_pct_d90_bias_pct <- bsBiasFun("km_naive_pct","bias_pct_90","naive_bs_pct")
bs_wtd_d90_bias_pct <- bsBiasFun("km_wtd","bias_pct_90","bs_wtd")

wrFun(naive_sample_d90_bias_pct, "/naive_sample_d90_bias_pct")
wrFun(bs_naive_pct_d90_bias_pct, "/bs_naive_pct_d90_bias_pct")
wrFun(bs_wtd_d90_bias_pct, "/bs_wtd_d90_bias_pct")

# d365
d365_obs <- sample_km_obs$sample_times_obs$d365$surv
naive_sample_d365_bias <- biasFun(sample_km_pred_naive,"sample_times_pred","d365",d365_obs)
wrFun(naive_sample_d365_bias, "/naive_sample_d365_bias")

bs_naive_pct_d365_bias <- bsBiasFun("km_naive_pct","bias_365","naive_bs_pct")
wrFun(bs_naive_pct_d365_bias, "/bs_naive_pct_d365_bias")
bs_wtd_d365_bias <- bsBiasFun("km_wtd","bias_365","bs_wtd")
wrFun(bs_wtd_d365_bias, "/bs_wtd_d365_bias")
# d365 pct bias
naive_sample_d365_bias_pct <- pctBiasFun(sample_km_pred_naive,"sample_times_pred","d365",d365_obs)
bs_naive_pct_d365_bias_pct <- bsBiasFun("km_naive_pct","bias_pct_365","naive_bs_pct")
bs_wtd_d365_bias_pct <- bsBiasFun("km_wtd","bias_pct_365","bs_wtd")
wrFun(naive_sample_d365_bias_pct, "/naive_sample_d365_bias_pct")
wrFun(bs_naive_pct_d365_bias_pct, "/bs_naive_pct_d365_bias_pct")
wrFun(bs_wtd_d365_bias_pct, "/bs_wtd_d365_bias_pct")

# median
median_obs <- as.data.frame.matrix(t(sample_km_obs$sample_median_obs$quantile))
naive_sample_med_bias <- biasFun(sample_km_pred_naive,"sample_median_pred","quantile", median_obs)
naive_sample_med_bias <- naive_sample_med_bias %>% mutate_at(vars(-alg),funs(round(.,1)))

bs_naive_pct_med_bias <- bsBiasFun("km_naive_pct","bias_med","naive_bs_pct")
bs_wtd_med_bias <- bsBiasFun("km_wtd","bias_med","bs_wtd")

wrFun(naive_sample_med_bias, "/naive_sample_med_bias")
wrFun(bs_naive_pct_med_bias, "/bs_naive_pct_med_bias")
wrFun(bs_wtd_med_bias, "/bs_wtd_med_bias")

# median pct bias
naive_sample_med_bias_pct <- pctBiasFun(sample_km_pred_naive,"sample_median_pred","quantile", median_obs)
bs_naive_pct_med_bias_pct <- bsBiasFun("km_naive_pct","bias_pct_med","naive_bs_pct")
bs_wtd_med_bias_pct <- bsBiasFun("km_wtd","bias_pct_med","bs_wtd")
wrFun(naive_sample_med_bias_pct, "/naive_sample_med_bias_pct")
wrFun(bs_naive_pct_med_bias_pct, "/bs_naive_pct_med_bias_pct")
wrFun(bs_wtd_med_bias_pct, "/bs_wtd_med_bias_pct")

