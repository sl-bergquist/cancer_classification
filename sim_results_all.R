##################################################
# CancerCLAS Multiclass Prediction
# Simulation Results
# Table output
##################################################
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
library(xtable)

base_path <- "~/base_path"
save_path <- "~/save_path"

sets <- c("set_1", "set_2", "set_3") 
pct_sets <- c("pct_set_1", "pct_set_2", "pct_set_3")
# 1 is accurate and certain
# 2 is accurate and uncertain
# 3 is inaccurate and uncertain

readFun <- function(file_name,set){
  files <- file.path(base_path, c(paste0(file_name, set, ".RDS")))
  list <- lapply(files, readRDS)
  return(list)
}

wrFun <- function(df, name){
  write.table(df, file=paste0(save_path, name, ".txt"),sep=",", quote=F, row.names=F)
}
boundsFun <- function(data){
  tmp <- data %>% select(setting,starts_with("lb"), starts_with("ub"))
  tmp <- tmp %>% mutate(mean_g1=paste0("(",lb_g1, ", ",ub_g1,")"),
                        mean_g2=paste0("(",lb_g2, ", ",ub_g2,")"),
                        mean_g3=paste0("(",lb_g3, ", ",ub_g3,")")) %>%
    select(setting, mean_g1, mean_g2, mean_g3)
  out <- data %>% select(setting, mean_g1, mean_g2, mean_g3) %>% rbind(., tmp) 
  out <- out %>% mutate(order=str_sub(setting,-1L)) %>% arrange(order, setting) %>%
    select(-order)
  return(out)
}
# coverage, threshdolds, ambiguity
cleanFun <- function(list,group,digits, colnames){
  tmp <- lapply(list, function(x)do.call(rbind.data.frame,x))
  tmp <- lapply(lapply(tmp, colMeans), round, digits=digits)
  tmp <- do.call(rbind.data.frame, tmp) # bind settings together in a df
  colnames(tmp) <- colnames # keep colnames
  tmp <- tmp %>% mutate(setting=paste0(group,sets)) %>% select(setting, everything())
  return(tmp)
}
colnames3 <- c("g1", "g2", "g3")
## coverage 
naive_cov <- readFun("simrep_class_coverage_naive_", sets)
colnames1 <- colnames(naive_cov[[1]][[1]])
naive_cov <- cleanFun(naive_cov, "naive_",2, colnames=colnames1)
wb_cov <- readFun("simrep_class_coverage_wtd_", sets)
wb_cov <- cleanFun(wb_cov, "wtd_",2, colnames=colnames1)
cov_df <- rbind.data.frame(naive_cov[1,], wb_cov[1,],naive_cov[2,], wb_cov[2,],
                           naive_cov[3,], wb_cov[3,])
wrFun(cov_df, "/coverage"); rm(naive_cov, wb_cov)
cov_df <- boundsFun(cov_df)
# Table B3: Bootstrap-based coverage
print(xtable(cov_df),include.rownames=F); rm(cov_df)

## amibiguity
wb_amb <- readFun("simrep_ambiguity_", sets)
ambFun <- function(list, measure){
  tmp <- lapply(list, function(x){lapply(x, function(y){y[[measure]]})})
  # turn each sim setting into a df of all simreps
  tmp <- lapply(tmp, function(x)do.call(rbind.data.frame,x))
  # get colmeans
  tmp_mean <- lapply(lapply(tmp, colMeans), round)
  # bind settings togther in df
  tmp_mean <- do.call(rbind.data.frame, tmp_mean)
  c_names <- c("l0", "lb_l0","ub_l0","l1","lb_l1","ub_l1","l2","lb_l2","ub_l2","l3","lb_l3","ub_l3")
  tmp_out <- setNames(tmp_mean, c_names)
  tmp_out <- tmp_out %>% mutate(setting=c("1","2","3")) %>% select(setting, everything())
  return(tmp_out)
}
wb_amb_all <- ambFun(wb_amb, "amb_all")
wb_amb_g1 <- ambFun(wb_amb, "amb_g1")
wb_amb_g2 <- ambFun(wb_amb, "amb_g2")
wb_amb_g3 <- ambFun(wb_amb, "amb_g3")
wrFun(wb_amb_all, "/amb_all"); wrFun(wb_amb_g1, "/amb_g1"); wrFun(wb_amb_g2, "/amb_g2"); wrFun(wb_amb_g3, "/amb_g3")
# ambiguity proportions 
# rather than doing within each simrep and then averaging proportions, do it across
  # rows based on the means across simulations
# proportions
# all
tmp_amb_all_pct <- wb_amb_all %>% mutate_at(vars(-setting), funs(./1000)) %>% mutate_at(vars(-setting), ~round(.,2))
boundsFun2 <- function(data){
  tmp <- data %>% select(setting,starts_with("lb"), starts_with("ub"))
  tmp <- tmp %>% mutate(l0=paste0("(",lb_l0, ", ",ub_l0,")"),
                        l1=paste0("(",lb_l1, ", ",ub_l1,")"),
                        l2=paste0("(",lb_l2, ", ",ub_l2,")"),
                        l3=paste0("(",lb_l3, ", ",ub_l3,")"))%>%
    select(setting, l0, l1, l2, l3)
  out <- data %>% select(setting,l0, l1, l2, l3) %>% rbind(., tmp) 
  out <- out %>% mutate(order=str_sub(setting,-1L)) %>% arrange(order, setting) %>%
    select(-order)
  return(out)
}
tmp_amb_all_pct <- boundsFun2(tmp_amb_all_pct)
print(xtable(tmp_amb_all_pct), include.rownames=F)
# by class
n_g1 <- (wb_amb_g1$l0 + wb_amb_g1$l1 + wb_amb_g1$l2 + wb_amb_g1$l3)[1]
n_g2 <- (wb_amb_g2$l0 + wb_amb_g2$l1 + wb_amb_g2$l2 + wb_amb_g2$l3)[1]
n_g3 <- (wb_amb_g3$l0 + wb_amb_g3$l1 + wb_amb_g3$l2 + wb_amb_g3$l3)[1]

acFun <- function(data_amb){
  tmp <- data %>% mutate(class=c(1,2,3))
  tmp_g1 <- tmp %>% filter(class==1) %>% mutate_at(vars(-class), funs(./n_g1))
  tmp_g2 <- tmp %>% filter(class==2) %>% mutate_at(vars(-class), funs(./n_g2))
  tmp_g3 <- tmp %>% filter(class==3) %>% mutate_at(vars(-class), funs(./n_g3))
  tmp <- rbind.data.frame(tmp_g1, tmp_g2, tmp_g3)
  tmp <- tmp %>% select(class, everything())
  tmp <- round(tmp,2)
  return(tmp)
}

tmp_amb_g1_pct <- wb_amb_g1 %>% mutate_at(vars(-setting), funs(./n_g1)) %>% mutate_at(vars(-setting), ~round(.,2))
tmp_amb_g1_pct <- boundsFun2(tmp_amb_g1_pct)
print(xtable(tmp_amb_g1_pct),include.rownames=F)

tmp_amb_g2_pct <- wb_amb_g2 %>% mutate_at(vars(-setting), funs(./n_g2)) %>% mutate_at(vars(-setting), ~round(.,2))
tmp_amb_g2_pct <- boundsFun2(tmp_amb_g2_pct)
print(xtable(tmp_amb_g2_pct),include.rownames=F)

tmp_amb_g3_pct <- wb_amb_g3 %>% mutate_at(vars(-setting), funs(./n_g3)) %>% mutate_at(vars(-setting), ~round(.,2))
tmp_amb_g3_pct <- boundsFun2(tmp_amb_g3_pct)
print(xtable(tmp_amb_g3_pct),include.rownames=F)

## thresholds
wb_thlds <- readFun("simrep_thresholds_", sets)
wb_thlds <- cleanFun(wb_thlds, "wtd_",3, colnames=colnames3)
wrFun(wb_thlds,"/thresholds")
# Table B4: Weighted labeling thresholds.
print(xtable(wb_thlds, digits=3), include.rownames=F)

## class measures 
# naive sample
naive_sample <- readFun("simrep_sample_perf_naive_", sets)
cmSampleFun <- function(list, measure,group){
  m <- lapply(list, function(x){lapply(x, function(y){y[[measure]]})})
  # turn each sim setting into df of all simreps
  m <- lapply(m, function(x)do.call(rbind.data.frame, x))
  # get colmeans
  m_mean <- lapply(lapply(m, colMeans, na.rm=T), round, digits=2)
  # bind settings together in df
  m_mean <- do.call(rbind.data.frame, m_mean) 
  if(measure %in% c("sens","spec","ppv","acc")){
  c_tmp <- c("all", "g1", "g2", "g3")
  m_out <- setNames(m_mean, c_tmp)
  # add in lb and ub cols, which will be NA -- in order to emphasize no-bs method
  m_out <- m_out %>% mutate(setting=paste0(group,sets),
                            lb_all=NA, ub_all=NA, lb_g1=NA, ub_g1=NA,
                            lb_g2=NA, ub_g2=NA, lb_g3=NA, ub_g3=NA) %>% 
    select(setting, all, lb_all, ub_all, g1, lb_g1, ub_g1, g2, lb_g2, ub_g2, 
           g3, lb_g3, ub_g3)
  } else {
    c_tmp <- c("g1","g2","g3")
    m_out <- setNames(m_mean, c_tmp)
    m_out <- m_out %>% mutate(setting=paste0(group,sets),
                              lb_g1=NA, ub_g1=NA,lb_g2=NA, ub_g2=NA, 
                              lb_g3=NA, ub_g3=NA) %>% 
      select(setting, g1, lb_g1, ub_g1, g2, lb_g2, ub_g2, g3, lb_g3, ub_g3)
    }
  
  return(m_out)
}
cm_list <- c("sens", "spec","ppv","acc")
naive_sample_cm <- lapply(cm_list,cmSampleFun, list=naive_sample, group="naive_")
names(naive_sample_cm) <- cm_list
# bs versions
naive_pct_cm <- readFun("simrep_class_measures_naive_", pct_sets)
wb_cm <- readFun("simrep_class_measures_wtd_", sets)
cmFun <- function(list,measure,group){
  # getting to specific measures at 3rd level of nesting (avoiding purrr b/c breaks)
  # get just first row -- naive version produces 500 rows and haven't fixed it yet
  m <- lapply(list, function(x){lapply(x, function(y){y[[measure]]})})
  # turn each sim setting into df of all simreps
  m <- lapply(m, function(x)do.call(rbind.data.frame, x))
  # get colmeans
  m <- lapply(lapply(m, colMeans, na.rm=T),round,digits=2)
  m <- do.call(rbind.data.frame, m) 
  
  if(measure %in% c("sens","spec","ppv","acc")){
    colnames <- c("all","lb_all","ub_all","g1", "lb_g1", "ub_g1","g2","lb_g2","ub_g2","g3","lb_g3","ub_g3")
  } else{
    colnames <- c("g1", "lb_g1", "ub_g1","g2","lb_g2","ub_g2","g3","lb_g3","ub_g3")
  }
  colnames(m) <- colnames
  # bind settings together in df and keep colnames  
  m <- m %>% mutate(setting=paste0(group,sets)) %>% select(setting, everything())
  return(m)
}
naive_pct_cm <- lapply(cm_list, cmFun, list=naive_pct_cm, group="naive_xbs_")
names(naive_pct_cm) <- cm_list
wb_cm <- lapply(cm_list, cmFun, list=wb_cm,group="wtd_bs_")
names(wb_cm) <- cm_list

## combine into df
sens_df <- rbind.data.frame(naive_sample_cm$sens, naive_pct_cm$sens, wb_cm$sens)
spec_df <- rbind.data.frame(naive_sample_cm$spec, naive_pct_cm$spec, wb_cm$spec)
ppv_df <- rbind.data.frame(naive_sample_cm$ppv, naive_pct_cm$ppv, wb_cm$ppv)
acc_df <- rbind.data.frame(naive_sample_cm$acc, naive_pct_cm$acc, wb_cm$acc)
wrFun(sens_df, "/sens"); wrFun(spec_df, "/spec"); wrFun(ppv_df, "/ppv"); wrFun(acc_df, "/acc")

boundsFun3 <- function(data){
  tmp <- data %>% select(setting,starts_with("lb"), starts_with("ub"))
  tmp <- tmp %>% mutate(all=paste0("(",lb_all, ", ",ub_all,")"),
                        g1=paste0("(",lb_g1, ", ",ub_g1,")"),
                        g2=paste0("(",lb_g2, ", ",ub_g2,")"),
                        g3=paste0("(",lb_g3, ", ",ub_g3,")"))%>%
    select(setting, all, g1, g2, g3)
  out <- data %>% select(setting,all, g1, g2, g3) %>% rbind(., tmp) 
  out <- out %>% mutate(order=str_sub(setting,-1L)) %>% arrange(order, setting) %>%
    select(-order)
  return(out)
}
sens_df <- boundsFun3(sens_df)
spec_df <- boundsFun3(spec_df)
ppv_df <- boundsFun3(ppv_df)
acc_df <- boundsFun3(acc_df)

# Table B5: Classification accuracy
print(xtable(acc_df), include.rownames=F)
# Table B6: Classification sensitivity
print(xtable(sens_df), include.rownames=F)
# Table B7: Classification specificity
print(xtable(spec_df), include.rownames=F)
# Table B8: Classification PPV
print(xtable(ppv_df), include.rownames=F)

## class counts
# naive sample
naive_sample <- readFun("simrep_sample_perf_naive_", sets)
cc_list <- c("tp", "tn", "fp", "fn")
naive_sample_cc <- lapply(cc_list,cmSampleFun, list=naive_sample, group="naive_")
names(naive_sample_cc) <- cc_list
# bs versions
naive_pct_cc <- readFun("simrep_class_measures_naive_", pct_sets)
wb_cc <- readFun("simrep_class_measures_wtd_", sets)
naive_pct_cc <- lapply(cc_list, cmFun, list=naive_pct_cc, group="naive_xbs_")
names(naive_pct_cc) <- cc_list
wb_cc <- lapply(cc_list, cmFun, list=wb_cc,group="wtd_bs_")
names(wb_cc) <- cc_list
## write each of these to text -- make as dfs and combine
tp_df <- rbind.data.frame(naive_sample_cc$tp,naive_pct_cc$tp, wb_cc$tp)
tn_df <- rbind.data.frame(naive_sample_cc$tn,naive_pct_cc$tn, wb_cc$tn)
fp_df <- rbind.data.frame(naive_sample_cc$fp,naive_pct_cc$fp, wb_cc$fp)
fn_df <- rbind.data.frame(naive_sample_cc$fn,naive_pct_cc$fn, wb_cc$fn)
 wrFun(tp_df, "/tp"); wrFun(tn_df, "/tn"); wrFun(fp_df, "/fp"); wrFun(fn_df, "/fn")
boundsFun4 <- function(data){
      tmp <- data %>% select(setting,starts_with("lb"), starts_with("ub"))
    tmp <- tmp %>% mutate(g1=paste0("(",lb_g1, ", ",ub_g1,")"),
                          g2=paste0("(",lb_g2, ", ",ub_g2,")"),
                          g3=paste0("(",lb_g3, ", ",ub_g3,")")) %>%
      select(setting, g1, g2, g3)
    out <- data %>% select(setting, g1, g2, g3) %>% rbind(., tmp) 
    out <- out %>% mutate(order=str_sub(setting,-1L)) %>% arrange(order, setting) %>%
      select(-order)
    return(out)
}
tp_df <- boundsFun4(tp_df) 
tn_df <- boundsFun4(tn_df)
fp_df <- boundsFun4(fp_df)
fn_df <- boundsFun4(fn_df)
# Table B9: True positive counts
print(xtable(tp_df), include.rownames=F)
# Table B10: True negative counts
print(xtable(tn_df), include.rownames=F)
# Table B11: False positive counts
print(xtable(fp_df), include.rownames=F)
# Table B12: False negative counts
print(xtable(fn_df), include.rownames=F)


## survival -- obs
km_obs <- readFun("simrep_km_obs_", sets)
kmTimeFun <- function(list,group, sample,days){
  # getting to specific measures at 3rd level of nesting (avoiding purrr b/c breaks)
  times <- lapply(list, function(x){lapply(x, function(y){y[[paste0("sample_times_",sample)]]})})
  # separate out specific times
  t <- lapply(times, function(x){lapply(x, function(y){y[[days]]})})
  # turn each sim setting into df of all simreps
  surv <- lapply(t, function(x){lapply(x, function(y){y[["surv"]]})})
  surv <- lapply(surv, function(x)do.call(rbind.data.frame, x)) 
  lb <- lapply(t, function(x){lapply(x, function(y){y[["lb"]]})})
  lb <- lapply(lb, function(x)do.call(rbind.data.frame, x))
  ub <- lapply(t, function(x){lapply(x, function(y){y[["ub"]]})})
  ub <- lapply(ub, function(x)do.call(rbind.data.frame, x))

  # get colmeans
  surv_mean <- lapply(lapply(surv, colMeans, na.rm=T), round, digits=2)
  lb_mean <- lapply(lapply(lb, colMeans, na.rm=T), round, digits=2)
  ub_mean <- lapply(lapply(ub, colMeans, na.rm=T), round, digits=2)
  # sd of mean survival and bounds across the simreps
  surv_sd <- lapply(lapply(surv, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
  lb_sd <- lapply(lapply(lb, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
  ub_sd <- lapply(lapply(ub, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
  
  # bind surv and surv sd together
  surv <- do.call(rbind.data.frame, surv_mean); surv_sd <- do.call(rbind.data.frame, surv_sd)
  lb <- do.call(rbind.data.frame, lb_mean); lb_sd <- do.call(rbind.data.frame, lb_sd)
  ub <- do.call(rbind.data.frame, ub_mean); ub_sd <- do.call(rbind.data.frame, ub_sd)
  
  # bind bounds and sds together
  colnames <- c("g1","g2","g3")
  colnames(surv) <- colnames; colnames(lb) <- paste0("lb_",colnames); colnames(ub) <- paste0("ub_",colnames)
  colnames_sd <- c("g1_sd","g2_sd","g3_sd")
  colnames(surv_sd) <- colnames_sd; colnames(lb_sd) <- paste0("lb_", colnames_sd); colnames(ub_sd) <- paste0("ub_",colnames_sd)
  surv <- cbind.data.frame(surv, surv_sd)
  bounds <- cbind.data.frame(lb, lb_sd, ub, ub_sd)
  

  surv <- surv %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  bounds  <- bounds %>% mutate(setting=paste0(group,sets)) %>% select(setting, lb_g1, lb_g1_sd, lb_g2, lb_g2_sd, lb_g3, lb_g3_sd,
                                                                      ub_g1, ub_g1_sd, ub_g2, ub_g2_sd, ub_g3, ub_g3_sd)
  out <- list(surv=surv, bounds=bounds)
  return(out)
}
km_obs_d30 <- kmTimeFun(km_obs, "obs_","obs", "d30")
km_obs_d90 <- kmTimeFun(km_obs, "obs_","obs", "d90")
km_obs_d365 <- kmTimeFun(km_obs, "obs_","obs", "d365")
kmMedFun <- function(list, group, sample,sets){
   median <- lapply(list, function(x){lapply(x, function(y){y[[paste0("sample_median_",sample)]]})})
  # what do I want from the median?
   # match structure of d30 etc: list of two (mean, bounds)
   surv <- lapply(median, function(x){lapply(x, function(y){y[["quantile"]]})})
   surv <- lapply(surv, function(x){lapply(x, t)})
   surv <- lapply(surv, function(x)do.call(rbind.data.frame, x)) 
   lb <- lapply(median, function(x){lapply(x, function(y){y[["lower"]]})})
   lb <- lapply(lb, function(x){lapply(x, t)})
   lb <- lapply(lb, function(x)do.call(rbind.data.frame, x)) 
   ub <- lapply(median, function(x){lapply(x, function(y){y[["upper"]]})})
   ub <- lapply(ub, function(x){lapply(x, t)})
   ub <- lapply(ub, function(x)do.call(rbind.data.frame, x)) 
   
   # get colmeans
   surv_mean <- lapply(lapply(surv, colMeans, na.rm=T), round, digits=2)
   lb_mean <- lapply(lapply(lb, colMeans, na.rm=T), round, digits=2)
   ub_mean <- lapply(lapply(ub, colMeans, na.rm=T), round, digits=2)
   # sd of mean survival and bounds across the simreps
   surv_sd <- lapply(lapply(surv, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
   lb_sd <- lapply(lapply(lb, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
   ub_sd <- lapply(lapply(ub, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
   
   # bind surv and surv sd together
   surv <- do.call(rbind.data.frame, surv_mean); surv_sd <- do.call(rbind.data.frame, surv_sd)
   lb <- do.call(rbind.data.frame, lb_mean); lb_sd <- do.call(rbind.data.frame, lb_sd)
   ub <- do.call(rbind.data.frame, ub_mean); ub_sd <- do.call(rbind.data.frame, ub_sd)
   
   # bind bounds and sds together
   colnames <- c("g1","g2","g3")
   colnames(surv) <- colnames; colnames(lb) <- paste0("lb_",colnames); colnames(ub) <- paste0("ub_",colnames)
   colnames_sd <- c("g1_sd","g2_sd","g3_sd")
   colnames(surv_sd) <- colnames_sd; colnames(lb_sd) <- paste0("lb_", colnames_sd); colnames(ub_sd) <- paste0("ub_",colnames_sd)
   surv <- cbind.data.frame(surv, surv_sd)
   bounds <- cbind.data.frame(lb, lb_sd, ub, ub_sd)
   
   surv <- surv %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
   bounds  <- bounds %>% mutate(setting=paste0(group,sets)) %>% select(setting, lb_g1, lb_g1_sd, lb_g2, lb_g2_sd, lb_g3, lb_g3_sd,
                                                                       ub_g1, ub_g1_sd, ub_g2, ub_g2_sd, ub_g3, ub_g3_sd)
   out <- list(surv=surv, bounds=bounds)
   return(out)
   
}
km_obs_median <- kmMedFun(km_obs,"obs_","obs", sets=sets)
# output: mean, lb, ub
formatFun <- function(data){
  out <- cbind(data$surv[,c("setting","g1")],
               data$bounds[,c("lb_g1","ub_g1")],
               data$surv[,c("g2")], data$bounds[,c("lb_g2","ub_g2")],
               data$surv[,c("g3")], data$bounds[,c("lb_g3","ub_g3")])
  colnames(out) <- c("setting", "g1", "lb_g1","ub_g1","g2","lb_g2","ub_g2","g3","lb_g3","ub_g3")
  return(out)
}
# median
km_obs_median_out <- formatFun(km_obs_median)
wrFun(km_obs_median_out,"/median_obs")
# d30
km_obs_d30_out <- formatFun(km_obs_d30)
wrFun(km_obs_d30_out,"/d30_obs")
# d90
km_obs_d90_out <- formatFun(km_obs_d90)
wrFun(km_obs_d90_out,"/d90_obs")
# d9365
km_obs_d365_out <- formatFun(km_obs_d365)
wrFun(km_obs_d365_out,"/d365_obs")

# survival -- predicted
# naive sample, no bs
naive_sample <- readFun("simrep_km_val_pred_sample_", sets)
# new function b/c here have to deal with setting 3 missing class 3
# honestly should just focus the below fun on setting 3
kmSample3Fun <- function(list,sample){
  median <- lapply(list, function(x){lapply(x, function(y){y[[paste0("sample_median_",sample)]]})})
  median <- lapply(median, function(x){lapply(x, t)})
  median <- median[[3]]
  # just doing median first but then can add in time points to clean in one go
  times <- lapply(list, function(x){lapply(x, function(y){y[[paste0("sample_times_",sample)]]})})
  times <- times[[3]]
  d30 <- lapply(times, `[[`,"d30")
  d90 <- lapply(times, `[[`,"d90")
  d365 <- lapply(times, `[[`,"d365")
  
  # check simreps -- see if there are fewer than 3
  # this is lazy and doens't distinguish btw classes
  srchkFun <- function(data,measure){
    tmp <- do.call(rbind.data.frame,lapply(lapply(data,  `[[`,measure),length))
    colnames(tmp) <- "count"
    out <- as.list(which(tmp$count<3))
    return(out)
  }
  problem_reps_median <- srchkFun(median, 1)
  problem_reps_time <- srchkFun(d30, "surv") # pretty sure these are the same but whatever
  
  # median survival time 
  for(i in problem_reps_median){
    median[[i]][[1]][3] <- NA
    median[[i]][[2]][3] <- NA
    median[[i]][[3]][3] <- NA
  }
  # then make all reps of median just numerics
  median <- lapply(median, function(x){lapply(x, as.numeric)})
  median <- lapply(median, setNames, c("quantile","lower","upper"))
  
  for(i in problem_reps_time){
    d30[[i]]$surv[3] <- NA
    d30[[i]]$lb[3] <- NA
    d30[[i]]$ub[3] <- NA
    d90[[i]]$surv[3] <- NA
    d90[[i]]$lb[3] <- NA
    d90[[i]]$ub[3] <- NA
    d365[[i]]$surv[3] <- NA
    d365[[i]]$lb[3] <- NA
    d365[[i]]$ub[3] <- NA
  }

# then just output to replace the third setting and can use other funs as before?
  sample_times_pred <- vector("list",1000)
  for(i in 1:1000){
    sample_times_pred[[i]]$d30 <- d30[[i]]
    sample_times_pred[[i]]$d90 <- d90[[i]]
    sample_times_pred[[i]]$d365 <- d365[[i]]
  }

  out <- list(sample_times_pred=sample_times_pred, sample_median_pred=median)
  
  return(out)
}
tmp <- kmSample3Fun(naive_sample, "pred")
tmp_list <- vector("list",1000)
for(i in 1:1000){
  tmp_list[[i]]$sample_times_pred <- tmp$sample_times_pred[[i]]
  tmp_list[[i]]$sample_median_pred <- tmp$sample_median_pred[[i]]
}
for(i in 1:1000){
  naive_sample[[3]][[i]] <- tmp_list[[i]]
}

naive_sample_km_d30 <- kmTimeFun(naive_sample, "naive_pred_","pred", "d30")
naive_sample_km_d90 <- kmTimeFun(naive_sample, "naive_pred_","pred", "d90")
naive_sample_km_d365 <- kmTimeFun(naive_sample, "naive_pred_","pred", "d365")
# time structure works fine, it's just median being a pain
# takes care of sets 1 and 2 fine
tmp_12 <- kmMedFun(naive_sample[1:2], "naive_pred_","pred", sets=sets[1:2])
# function to just deal with the third setting
kmMedSample3Fun <- function(list, sample){
  median <- lapply(list, function(x){lapply(x, function(y){y[[paste0("sample_median_",sample)]]})})
  median <- median[[3]] # this function is just for the third scenario
  # match structure of d30 etc: list of two (mean, bounds)
  surv <- lapply(median, `[[`,1) 
  surv <- do.call(rbind.data.frame,surv)
  lb <- lapply(median, `[[`,2)
  lb <- do.call(rbind.data.frame, lb)
  ub <- lapply(median, `[[`,3)
  ub <- do.call(rbind.data.frame, ub)
  
  # get colmeans
  surv_mean <- round(colMeans(surv, na.rm=T), digits=2)
  lb_mean <- round(colMeans(surv, na.rm=T),digits=2)
  ub_mean <- round(colMeans(surv, na.rm=T),digits=2)
  # sd of mean survival and bounds across the simreps
  surv_sd <- round(apply(surv,2,sd, na.rm=T), digits=2)
  lb_sd <- round(apply(lb, 2, sd, na.rm=T),digits=2)
  ub_sd <- round(apply(ub, 2, sd, na.rm=T),digits=2)
  
  surv <- as.vector(c(surv_mean[1],surv_sd[1],surv_mean[2],surv_sd[2],surv_mean[3],surv_sd[3]))
  surv <- rbind.data.frame(surv)
  c_names <- c("g1","g1_sd", "g2","g2_sd","g3","g3_sd")
  colnames(surv) <- c_names
  surv <- surv %>% mutate(setting=3) %>% select(setting, everything())
  lb <- as.vector(c(lb_mean[1],lb_sd[1],lb_mean[2], lb_sd[2], lb_mean[3], lb_sd[3]))
  lb <- rbind.data.frame(lb)
  colnames(lb) <- paste0("lb_",c_names)
  ub <- as.vector(c(ub_mean[1],ub_sd[1],ub_mean[2],ub_sd[2],ub_mean[3],ub_sd[3]))
  ub <- rbind.data.frame(ub)
  colnames(ub) <- paste0("ub_", c_names)
  bounds <- cbind.data.frame(lb, ub)
  bounds <- bounds %>% mutate(setting=3) %>% select(setting, everything())
  
  out <- list(surv=surv, bounds=bounds)
  return(out)
}
tmp_3 <- kmMedSample3Fun(naive_sample,"pred")
# combine settings 1 and 2 with setting 3 
naive_sample_km_median <- NULL
naive_sample_km_median$surv <- rbind.data.frame(tmp_12$surv, tmp_3$surv)
naive_sample_km_median$bounds <- rbind.data.frame(tmp_12$bounds, tmp_3$bounds)
rm(tmp_12, tmp_3)

# bs versions -- need new funs b/c structure is different
naive_pct_km <- readFun("simrep_km_naive_", pct_sets)
wb_km <- readFun("simrep_km_wtd_", sets)
kmBSFun <- function(list, group, days){
  tmp <- lapply(list, function(x){lapply(x, function(y){y[[days]]})})
  # combine into df
  tmp <- lapply(tmp, function(x)do.call(rbind.data.frame, x))
  # colmeans + sd
  tmp_means <- lapply(lapply(tmp, colMeans,na.rm=T),round, digits=2)
  tmp_sd <- lapply(lapply(tmp, function(x)apply(x, 2, sd, na.rm=T)), round,digits=2)
  # put in same structure as sample vesion
  surv <- lapply(tmp_means, `[`,c("g1","g2","g3")); surv_sd <- lapply(tmp_sd, `[`, c("g1","g2","g3"))
  lb <- lapply(tmp_means, `[`,c("lb_g1","lb_g2","lb_g3")); lb_sd <- lapply(tmp_sd, `[`,c("lb_g1","lb_g2","lb_g3"))
  ub <- lapply(tmp_means, `[`,c("ub_g1","ub_g2","ub_g3")); ub_sd <- lapply(tmp_sd, `[`,c("ub_g1","ub_g2","ub_g3"))
  # bind surv and sd together
  surv <- do.call(rbind.data.frame, surv); surv_sd <- do.call(rbind.data.frame, surv_sd)
  lb <- do.call(rbind.data.frame, lb); lb_sd <- do.call(rbind.data.frame, lb_sd)
  ub <- do.call(rbind.data.frame, ub); ub_sd <- do.call(rbind.data.frame, ub_sd)
  
  # bind bounds and sds together
  colnames <- c("g1","g2","g3")
  colnames(surv) <- colnames; colnames(lb) <- paste0("lb_",colnames); colnames(ub) <- paste0("ub_",colnames)
  colnames_sd <- c("g1_sd","g2_sd","g3_sd")
  colnames(surv_sd) <- colnames_sd; colnames(lb_sd) <- paste0("lb_", colnames_sd); colnames(ub_sd) <- paste0("ub_",colnames_sd)
  surv <- cbind.data.frame(surv, surv_sd)
  bounds <- cbind.data.frame(lb, lb_sd, ub, ub_sd)
  
  
  surv <- surv %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  bounds  <- bounds %>% mutate(setting=paste0(group,sets)) %>% select(setting, lb_g1, lb_g1_sd, lb_g2, lb_g2_sd, lb_g3, lb_g3_sd,
                                                                      ub_g1, ub_g1_sd, ub_g2, ub_g2_sd, ub_g3, ub_g3_sd)
  out <- list(surv=surv, bounds=bounds)
  return(out)
  
}
days_list <- c("d30","d90","d365","median")
naive_pct_km_list <- lapply(days_list, kmBSFun, list=naive_pct_km, group="naive_xbs_")
names(naive_pct_km_list) <- days_list

wb_km_list <- lapply(days_list, kmBSFun, list=wb_km, group="wtd_bs_")
names(wb_km_list) <- days_list

# median
km_pred_median <- formatFun(naive_sample_km_median)
naive_pct_km_median <- formatFun(naive_pct_km_list$median)
wb_km_median <- formatFun(wb_km_list$median)
median_out <- rbind(km_pred_median[1,],naive_pct_km_median[1,],wb_km_median[1,],
                    km_pred_median[2,],naive_pct_km_median[2,],wb_km_median[2,],
                    km_pred_median[3,],naive_pct_km_median[3,],wb_km_median[3,])
wrFun(median_out, "/median_pred")
median_out <- boundsFun4(median_out)
print(xtable(median_out), include.rownames=F)
# d30
km_pred_d30 <- formatFun(naive_sample_km_d30)
naive_pct_km_d30 <- formatFun(naive_pct_km_list$d30)
wb_km_d30 <- formatFun(wb_km_list$d30)
d30_out <- rbind(km_pred_d30[1,],naive_pct_km_d30[1,],wb_km_d30[1,],
                    km_pred_d30[2,],naive_pct_km_d30[2,],wb_km_d30[2,],
                 km_pred_d30[3,],naive_pct_km_d30[3,],wb_km_d30[3,])
wrFun(d30_out, "/d30_pred")
d30_out <- boundsFun4(d30_out)
# d90
km_pred_d90 <- formatFun(naive_sample_km_d90)
naive_pct_km_d90 <- formatFun(naive_pct_km_list$d90)
wb_km_d90 <- formatFun(wb_km_list$d90)
d90_out <- rbind(km_pred_d90[1,],naive_pct_km_d90[1,],wb_km_d90[1,],
                 km_pred_d90[2,],naive_pct_km_d90[2,],wb_km_d90[2,],
                 km_pred_d90[3,],naive_pct_km_d90[3,],wb_km_d90[3,])
wrFun(d90_out, "/d90_pred")
d90_out <- boundsFun4(d90_out)
# d365
km_pred_d365 <- formatFun(naive_sample_km_d365)
naive_pct_km_d365 <- formatFun(naive_pct_km_list$d365)
wb_km_d365 <- formatFun(wb_km_list$d365)
d365_out <- rbind(km_pred_d365[1,],naive_pct_km_d365[1,],wb_km_d365[1,],
                 km_pred_d365[2,],naive_pct_km_d365[2,],wb_km_d365[2,],
                 km_pred_d365[3,],naive_pct_km_d365[3,],wb_km_d365[3,])
wrFun(d365_out, "/d365_pred")
d365_out <- boundsFun4(d365_out)


# survival: bias -- comparing predicted to obs
km_obs <- readFun("simrep_km_obs_", sets)
# use naive sample from above that has been treated with kmSample3Fun
kmTimeBiasFun <- function(list_obs, list,group, days){
  # getting to specific measures at 3rd level of nesting (avoiding purrr b/c breaks)
  cleanUpFun <- function(data, sample){
    times <- lapply(data, function(x){lapply(x, function(y){y[[paste0("sample_times_",sample)]]})})
    t <- lapply(times, function(x){lapply(x, function(y){y[[days]]})})
    
    surv <- lapply(t, function(x){lapply(x, function(y){y[["surv"]]})})
    surv <- lapply(surv, function(x)do.call(rbind.data.frame, x)) 
    lb <- lapply(t, function(x){lapply(x, function(y){y[["lb"]]})})
    lb <- lapply(lb, function(x)do.call(rbind.data.frame, x))
    ub <- lapply(t, function(x){lapply(x, function(y){y[["ub"]]})})
    ub <- lapply(ub, function(x)do.call(rbind.data.frame, x))
    colnames <- c("g1","g2","g3")
    surv <- lapply(surv, setNames, colnames); lb <- lapply(lb, setNames, colnames); ub <- lapply(ub, setNames, colnames)
    out <- list(surv=surv, lb=lb, ub=ub)
    return(out)
  }
  obs_out <- cleanUpFun(data=list_obs, sample="obs")
  pred_out <- cleanUpFun(data=list, sample="pred")
  
  # calculate bias: first doing means
  biasFun <- function(pred, obs, measure, set){
    tmp1 <- pred[[measure]][[set]]["g1"] - obs[[measure]][[set]]["g1"]
    tmp2 <- pred[[measure]][[set]]["g2"] - obs[[measure]][[set]]["g2"]
    tmp3 <- pred[[measure]][[set]]["g3"] - obs[[measure]][[set]]["g3"]
    out <- cbind(tmp1, tmp2, tmp3)
    colnames(out) <- c("g1","g2","g3")
    return(out)
  }
  surv_1 <- biasFun(pred=pred_out,obs=obs_out, measure="surv",set=1)
  surv_2 <- biasFun(pred=pred_out,obs=obs_out, measure="surv",set=2)
  surv_3 <- biasFun(pred=pred_out, obs=obs_out, measure="surv",set=3)
  lb_1 <- biasFun(pred=pred_out,obs=obs_out, measure="lb",set=1)
  lb_2 <- biasFun(pred=pred_out,obs=obs_out, measure="lb",set=2)
  lb_3 <- biasFun(pred=pred_out, obs=obs_out, measure="lb",set=3)
  ub_1 <- biasFun(pred=pred_out,obs=obs_out, measure="ub",set=1)
  ub_2 <- biasFun(pred=pred_out,obs=obs_out, measure="ub",set=2)
  ub_3 <- biasFun(pred=pred_out,obs=obs_out, measure="ub",set=3)
  
  surv_list <- list(surv_1, surv_2, surv_3); lb_list <- list(lb_1, lb_2, lb_3); ub_list <- list(ub_1,ub_2,ub_3)
  # get mean and sd bias
  surv_mean <- lapply(lapply(surv_list, colMeans, na.rm=T),round,digits=3)
  surv_sd <- lapply(lapply(surv_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  lb_mean <- lapply(lapply(lb_list, colMeans, na.rm=T),round,digits=3)
  lb_sd <- lapply(lapply(lb_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  ub_mean <- lapply(lapply(ub_list, colMeans, na.rm=T),round,digits=3)
  ub_sd <- lapply(lapply(ub_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  surv_lb <- lapply(surv_list,function(x)apply(x,2,quantile,.025,na.rm=T ))
  surv_ub <- lapply(surv_list,function(x)apply(x,2,quantile,.975,na.rm=T))
  
  # bind settings together
  surv_mean <- do.call(rbind.data.frame, surv_mean); surv_sd <- do.call(rbind.data.frame, surv_sd)
  surv_lb <- do.call(rbind.data.frame, surv_lb); surv_ub <- do.call(rbind.data.frame,surv_ub)
  lb_mean <- do.call(rbind.data.frame, lb_mean); lb_sd <- do.call(rbind.data.frame, lb_sd)
  ub_mean <- do.call(rbind.data.frame, ub_mean); ub_sd <- do.call(rbind.data.frame, ub_sd)
  
  c_tmp <- c("g1","g2","g3")
  sd_tmp <- c("g1_sd","g2_sd","g3_sd")
  colnames(surv_mean) <- c_tmp; colnames(lb_mean) <- c_tmp; colnames(ub_mean) <- c_tmp
  colnames(surv_sd) <- sd_tmp; colnames(lb_sd) <- sd_tmp; colnames(ub_sd) <- sd_tmp
  colnames(surv_lb) <- c_tmp; colnames(surv_ub) <- c_tmp
  surv_out <- cbind.data.frame(surv_mean, surv_sd); lb_out <- cbind.data.frame(lb_mean, lb_sd); ub_out <- cbind.data.frame(ub_mean, ub_sd)
  surv_out <- surv_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  lb_out <- lb_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  ub_out <- ub_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  surv_lb_out <- surv_lb %>% mutate(setting=paste0(group,sets)) %>% select(setting,g1,g2,g3)
  surv_ub_out <- surv_ub %>% mutate(setting=paste0(group, sets)) %>% select(setting,g1,g2,g3)
  out <- list(surv=surv_out,surv_lb=surv_lb_out,surv_ub=surv_ub_out, lb=lb_out, ub=ub_out)
  
  return(out)
}
naive_sample_bias_d30 <- kmTimeBiasFun(list_obs=km_obs, list=naive_sample,group="naive_",days="d30")  
naive_sample_bias_d90 <- kmTimeBiasFun(list_obs=km_obs, list=naive_sample,group="naive_",days="d90")  
naive_sample_bias_d365 <- kmTimeBiasFun(list_obs=km_obs, list=naive_sample,group="naive_",days="d365")  

# median
kmMedSampleBiasFun <- function(list_obs, list, group, sets){
  cleanUpFun <- function(data, sample){
    median <- lapply(data, function(x){lapply(x, function(y){y[[paste0("sample_median_",sample)]]})})
    # what do I want from the median?
    # match structure of d30 etc: list of two (mean, bounds)
    surv <- lapply(median, function(x){lapply(x, function(y){y[["quantile"]]})})
    surv <- lapply(surv, function(x){lapply(x, t)})
    surv <- lapply(surv, function(x)do.call(rbind.data.frame, x)) 
    lb <- lapply(median, function(x){lapply(x, function(y){y[["lower"]]})})
    lb <- lapply(lb, function(x){lapply(x, t)})
    lb <- lapply(lb, function(x)do.call(rbind.data.frame, x)) 
    ub <- lapply(median, function(x){lapply(x, function(y){y[["upper"]]})})
    ub <- lapply(ub, function(x){lapply(x, t)})
    ub <- lapply(ub, function(x)do.call(rbind.data.frame, x)) 
    colnames <- c("g1","g2","g3")
    surv <- lapply(surv, setNames, colnames); lb <- lapply(lb, setNames, colnames); ub <- lapply(ub, setNames, colnames)
    out <- list(surv=surv, lb=lb, ub=ub)
    return(out)
  }
  obs_out <- cleanUpFun(data=list_obs, sample="obs")
  pred_out <- cleanUpFun(data=list, sample="pred")
  # calculate bias
  biasFun <- function(pred, obs, measure, set){
    tmp1 <- pred[[measure]][[set]]["g1"] - obs[[measure]][[set]]["g1"]
    tmp2 <- pred[[measure]][[set]]["g2"] - obs[[measure]][[set]]["g2"]
    tmp3 <- pred[[measure]][[set]]["g3"] - obs[[measure]][[set]]["g3"]
    out <- cbind(tmp1, tmp2, tmp3)
    colnames(out) <- c("g1","g2","g3")
    return(out)
  }
  surv_1 <- biasFun(pred=pred_out,obs=obs_out, measure="surv",set=1)
  surv_2 <- biasFun(pred=pred_out,obs=obs_out, measure="surv",set=2)
  lb_1 <- biasFun(pred=pred_out,obs=obs_out, measure="lb",set=1)
  lb_2 <- biasFun(pred=pred_out,obs=obs_out, measure="lb",set=2)
  ub_1 <- biasFun(pred=pred_out,obs=obs_out, measure="ub",set=1)
  ub_2 <- biasFun(pred=pred_out,obs=obs_out, measure="ub",set=2)
  
  surv_list <- list(surv_1, surv_2); lb_list <- list(lb_1, lb_2); ub_list <- list(ub_1,ub_2)
  # get mean and sd bias
  surv_mean <- lapply(lapply(surv_list, colMeans, na.rm=T),round,digits=3)
  surv_sd <- lapply(lapply(surv_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  lb_mean <- lapply(lapply(lb_list, colMeans, na.rm=T),round,digits=3)
  lb_sd <- lapply(lapply(lb_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  ub_mean <- lapply(lapply(ub_list, colMeans, na.rm=T),round,digits=3)
  ub_sd <- lapply(lapply(ub_list, function(x)apply(x, 2, sd, na.rm=T)), round,digits=3)
  surv_lb <- lapply(surv_list,function(x)apply(x,2,quantile,.025,na.rm=T))
  surv_ub <- lapply(surv_list,function(x)apply(x,2,quantile,.975,na.rm=T))  
  
  # bind settings together
  surv_mean <- do.call(rbind.data.frame, surv_mean); surv_sd <- do.call(rbind.data.frame, surv_sd)
  lb_mean <- do.call(rbind.data.frame, lb_mean); lb_sd <- do.call(rbind.data.frame, lb_sd)
  ub_mean <- do.call(rbind.data.frame, ub_mean); ub_sd <- do.call(rbind.data.frame, ub_sd)
  surv_lb <- do.call(rbind.data.frame, surv_lb); surv_ub <- do.call(rbind.data.frame,surv_ub)
  
  c_tmp <- c("g1","g2","g3")
  sd_tmp <- c("g1_sd","g2_sd","g3_sd")
  colnames(surv_mean) <- c_tmp; colnames(lb_mean) <- c_tmp; colnames(ub_mean) <- c_tmp
  colnames(surv_sd) <- sd_tmp; colnames(lb_sd) <- sd_tmp; colnames(ub_sd) <- sd_tmp
  colnames(surv_ub) <- c_tmp; colnames(surv_lb) <- c_tmp
  surv_out <- cbind.data.frame(surv_mean, surv_sd); lb_out <- cbind.data.frame(lb_mean, lb_sd); ub_out <- cbind.data.frame(ub_mean, ub_sd)
  surv_out <- surv_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  lb_out <- lb_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  ub_out <- ub_out %>% mutate(setting=paste0(group,sets)) %>% select(setting, g1, g1_sd, g2, g2_sd, g3, g3_sd)
  surv_lb_out <- surv_lb %>% mutate(setting=paste0(group,sets)) %>% select(setting,g1,g2,g3)
  surv_ub_out <- surv_ub %>% mutate(setting=paste0(group, sets)) %>% select(setting,g1,g2,g3)
  out <- list(surv=surv_out,surv_lb=surv_lb_out,surv_ub=surv_ub_out, lb=lb_out, ub=ub_out)
  
  return(out)
}
tmp_12 <- kmMedSampleBiasFun(km_obs[1:2], naive_sample[1:2], group="naive_sample_", sets=sets[1:2])

kmMedSample3BiasFun <- function(list_obs, list, group){
  cleanUpFun <- function(data, sample){
  median <- lapply(data, function(x){lapply(x, function(y){y[[paste0("sample_median_",sample)]]})})
  median <- median[[3]]
  surv <- lapply(median, function(y){y[["quantile"]]})
  surv <- lapply(surv, function(x){lapply(x, t)})
  surv <- do.call(rbind.data.frame, surv) 
  lb <- lapply(median, function(y){y[["lower"]]})
  lb <- lapply(lb, function(x){lapply(x, t)})
  lb <- do.call(rbind.data.frame, lb)
  ub <- lapply(median, function(y){y[["upper"]]})
  ub <- lapply(ub, function(x){lapply(x, t)})
  ub <- do.call(rbind.data.frame, ub)
  c_names <- c("g1","g2","g3")
  colnames(surv) <- c_names; colnames(lb) <- c_names; colnames(ub) <- c_names
  out <- list(surv=surv, lb=lb, ub=ub)
  return(out)
  }
  obs_out <- cleanUpFun(data=list_obs, sample="obs")
  pred_out <- cleanUpFun(data=list, sample="pred")
  # calculate bias
  biasFun <- function(pred, obs, measure){
    tmp1 <- pred[[measure]]["g1"] - obs[[measure]]["g1"]
    tmp2 <- pred[[measure]]["g2"] - obs[[measure]]["g2"]
    tmp3 <- pred[[measure]]["g3"] - obs[[measure]]["g3"]
    out <- cbind(tmp1, tmp2, tmp3)
    colnames(out) <- c("g1","g2","g3")
    return(out)
  }
  surv_3 <- biasFun(pred=pred_out,obs=obs_out, measure="surv")
  lb_3 <- biasFun(pred=pred_out,obs=obs_out, measure="lb")
  ub_3 <- biasFun(pred=pred_out,obs=obs_out, measure="ub")

  # get mean and sd bias
  surv_mean <- round(colMeans(surv_3, na.rm=T),digits=3)
  surv_sd <- round(apply(surv_3,2,sd,na.rm=T),digits=3)
  surv_lb <- apply(surv_3,2,quantile,.025,na.rm=T )
  surv_ub <- apply(surv_3,2,quantile,.975,na.rm=T)
  
  lb_mean <- round(colMeans(lb_3,na.rm=T),digits=3)
  lb_sd <- round(apply(lb_3,2,sd,na.rm=T),digits=3)
  ub_mean <- round(colMeans(ub_3,na.rm=T),digits=3)
  ub_sd <- round(apply(ub_3,2,sd,na.rm=T),digits=3)

  # bind together
  c_names <- c("g1","g1_sd", "g2","g2_sd","g3","g3_sd")
  surv <- as.vector(c(surv_mean[1],surv_sd[1],surv_mean[2],surv_sd[2],surv_mean[3],surv_sd[3]))
  surv <- rbind.data.frame(surv)
  colnames(surv) <- c_names
  surv <- surv %>% mutate(setting=3) %>% select(setting, everything())
  
  lb <- as.vector(c(lb_mean[1],lb_sd[1],lb_mean[2],lb_sd[2],lb_mean[3],lb_sd[3]))
  lb <- rbind.data.frame(lb)
  colnames(lb) <- c_names
  lb <- lb %>% mutate(setting=3) %>% select(setting, everything())
  
  ub <- as.vector(c(ub_mean[1],ub_sd[1],ub_mean[2],ub_sd[2],ub_mean[3],ub_sd[3]))
  ub <- rbind.data.frame(ub)
  colnames(ub) <- c_names
  ub <- ub %>% mutate(setting=3) %>% select(setting, everything())
  
  surv_lb <- rbind.data.frame(surv_lb); surv_ub <- rbind.data.frame(surv_ub)
  colnames(surv_lb) <- c("g1","g2","g3"); colnames(surv_ub) <- c("g1","g2","g3")
  surv_lb <- surv_lb %>% mutate(setting=3) %>% select(setting, everything())
  surv_ub <- surv_ub %>% mutate(setting=3) %>% select(setting, everything())
  
  out <- list(surv=surv,surv_lb=surv_lb,surv_ub=surv_ub, lb=lb, ub=ub)
  
  return(out)
}
tmp_3 <- kmMedSample3BiasFun(list_obs=km_obs, list=naive_sample, group="naive_sample_")
# combine settings 1 and 2 with setting 3 
naive_sample_km_bias_median <- NULL
naive_sample_km_bias_median$surv <- rbind.data.frame(tmp_12$surv, tmp_3$surv)
naive_sample_km_bias_median$surv_lb <- rbind.data.frame(tmp_12$surv_lb, tmp_3$surv_lb)
naive_sample_km_bias_median$surv_ub <- rbind.data.frame(tmp_12$surv_ub, tmp_3$surv_ub)
naive_sample_km_bias_median$lb <- rbind.data.frame(tmp_12$lb, tmp_3$lb)
naive_sample_km_bias_median$ub <- rbind.data.frame(tmp_12$ub, tmp_3$ub)
rm(tmp_12, tmp_3)

# bs versions
naive_pct_km <- readFun("simrep_km_naive_", pct_sets)
wb_km <- readFun("simrep_km_wtd_", sets)
kmBiasBSFun <- function(list, group, measure){
  tmp <- lapply(list, function(x){lapply(x, function(y){y[[measure]]})})
  # combine into df
  tmp <- lapply(tmp, function(x)do.call(rbind.data.frame, x))
  pred_surv <- lapply(tmp,`[`,c("mean_g1","mean_g2","mean_g3"))
  pred_lb <- lapply(tmp,`[`,c("lb_g1","lb_g2","lb_g3")); pred_lb <- lapply(pred_lb, setNames, c("g1","g2","g3")) 
  pred_ub <- lapply(tmp,`[`,c("ub_g1","ub_g2","ub_g3")); pred_ub <- lapply(pred_ub, setNames, c("g1","g2","g3")) 

  # get colmeans -- replace INf with na first?
  pred_surv <- lapply(pred_surv, function(x){apply(x,2,function(y)ifelse(is.infinite(y),NA,y))})
  pred_lb <- lapply(pred_lb, function(x){apply(x,2,function(y)ifelse(is.infinite(y),NA,y))})
  pred_ub <- lapply(pred_ub, function(x){apply(x,2,function(y)ifelse(is.infinite(y),NA,y))})
  
  pred_surv <- do.call(rbind.data.frame,lapply(pred_surv, colMeans, na.rm=T))
  pred_lb <- do.call(rbind.data.frame, lapply(pred_lb, colMeans, na.rm=T))
  pred_ub <- do.call(rbind.data.frame, lapply(pred_ub, colMeans, na.rm=T))
  colnames(pred_surv) <- c("g1","g2","g3")
  colnames(pred_lb) <- c("g1","g2", "g3"); colnames(pred_ub) <- c("g1","g2","g3")
  
  # bind settings together in df and keep colnames
  pred_surv <- pred_surv %>% mutate(setting=paste0(group,sets)) %>% select(setting,everything())
  pred_lb <- pred_lb  %>% mutate(setting=paste0(group, sets)) %>% select(setting, everything())
  pred_ub <- pred_ub %>% mutate(setting=paste0(group, sets)) %>% select(setting, everything())
  
  out <- list(surv=pred_surv, surv_lb=pred_lb, surv_ub=pred_ub)
 
 return(out)
 
}
bias_list <- c("bias_30", "bias_90","bias_365","bias_med")
bias_pct_list <- c("bias_pct_30", "bias_pct_90","bias_pct_365","bias_pct_med")

# bias
naive_pct_km_bias_list <- lapply(bias_list, kmBiasBSFun, list=naive_pct_km,group="naive_bs_")
names(naive_pct_km_bias_list) <- bias_list
wb_km_bias_list <- lapply(bias_list, kmBiasBSFun, list=wb_km,group="wtd_bs_")
names(wb_km_bias_list) <- bias_list

outFun <- function(naive, measure1, measure2){
  naive <- naive[[measure1]] %>% select(setting, g1, g2, g3)
  out <- rbind.data.frame(naive[1,],
                                    naive_pct_km_bias_list[[measure2]][[measure1]][1,],
                                    wb_km_bias_list[[measure2]][[measure1]][1,],
                                    naive[2,],
                                    naive_pct_km_bias_list[[measure2]][[measure1]][2,],
                                    wb_km_bias_list[[measure2]][[measure1]][2,],
                                    naive[3,],
                                    naive_pct_km_bias_list[[measure2]][[measure1]][3,],
                                    wb_km_bias_list[[measure2]][[measure1]][3,])
  return(out)
}
# 30 days
d30_bias_mean <- outFun(naive_sample_bias_d30, "surv", "bias_30")
d30_bias_lb <- outFun(naive_sample_bias_d30,"surv_lb","bias_30")
d30_bias_ub <- outFun(naive_sample_bias_d30,"surv_ub","bias_30")
wrFun(d30_bias_mean, "/d30_bias_mean")
wrFun(d30_bias_lb, "/d30_bias_lb")
wrFun(d30_bias_ub, "/d30_bias_ub")
# 90 days
d90_bias_mean <- outFun(naive_sample_bias_d90, "surv", "bias_90")
d90_bias_lb <- outFun(naive_sample_bias_d90,"surv_lb","bias_90")
d90_bias_ub <- outFun(naive_sample_bias_d90,"surv_ub","bias_90")
wrFun(d90_bias_mean, "/d90_bias_mean")
wrFun(d90_bias_lb, "/d90_bias_lb")
wrFun(d90_bias_ub, "/d90_bias_ub")
# 365 days
d365_bias_mean <- outFun(naive_sample_bias_d365, "surv", "bias_365")
d365_bias_lb <- outFun(naive_sample_bias_d365,"surv_lb","bias_365")
d365_bias_ub <- outFun(naive_sample_bias_d365,"surv_ub","bias_365")
wrFun(d365_bias_mean, "/d365_bias_mean")
wrFun(d365_bias_lb, "/d365_bias_lb")
wrFun(d365_bias_ub, "/d365_bias_ub")
# median
med_bias_mean <- outFun(naive_sample_km_bias_median, "surv", "bias_med")
med_bias_lb <- outFun(naive_sample_km_bias_median,"surv_lb","bias_med")
med_bias_ub <- outFun(naive_sample_km_bias_median,"surv_ub","bias_med")
wrFun(med_bias_mean, "/median_bias_mean")
wrFun(med_bias_lb, "/median_bias_lb")
wrFun(med_bias_ub, "/median_bias_ub")
