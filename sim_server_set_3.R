#########################################################
# CancerCLAS Multiclass Prediction
# Simulate data from mvnorm
# Version for running on Server
# Naive No Variation, Naive Bootstrap, Weighted Bootstrap
#########################################################
options(scipen = 999)

library(MASS)
library(tidyverse)
library(glmnet)
library(nnet)
library(randomForest)
library(caret)
library(pROC)
library(survival)

# source misc_funs.R: contains function for calculating class measures
source("~/misc_funs_v3.R")

repname <- c("set_3") # inaccurate, uncertain
# Generate some example training data
set.seed(33)

#### Simulation Set Up ####
# set simulation repitions and sample sizes
n_dev <- 1000
n_val <- 1000
n_simrep <- 1000
# select number of bootstrap iterations
n_bs <- 500
# alpha for setting label thresholds
alpha <- .10
# objects to hold outcomes from simulations
simrep_thresholds <- vector(mode="list", n_simrep)
simrep_ambiguity <- vector(mode="list", n_simrep)
simrep_class_coverage_naive <- vector(mode="list", n_simrep)
simrep_class_coverage_wtd <- vector(mode="list", n_simrep)
simrep_class_measures_naive_pct <- vector(mode="list", n_simrep)
simrep_class_measures_wtd <- vector(mode="list", n_simrep)
simrep_km_obs <- vector(mode="list", n_simrep)
simrep_km_val_pred_sample <- vector(mode="list", n_simrep)
simrep_km_naive_pct <- vector(mode="list", n_simrep)
simrep_km_wtd <- vector(mode="list", n_simrep)
simrep_sample_perf_naive <- vector(mode="list", n_simrep)


#### Simrep Contents ####
system.time(
  for(r in 1:n_simrep){
    ## Draw Data ## 
    # features for prediction: uncorrelated
    X <- matrix(data=NA, nrow=n_dev+n_val, ncol=15)
    X[,1] <- rnorm(n_dev+n_val, 72,5) # age
    X[,2] <- rnorm(n_dev+n_val,45000,10000) # income
    X[,3] <- rnorm(n_dev+n_val,23,4) # Bmi-ish
    X[,4] <- rnorm(n_dev+n_val,70,5) # higher education
    X[,5] <- rnorm(n_dev+n_val,5,2) # normal 
    X[,6] <- rnorm(n_dev+n_val,0,1) # for survival
    X[,7] <- rbinom(n=n_dev+n_val,1,.5) # gender
    X[,8] <- rbinom(n=n_dev+n_val,1,.25) # surgery
    X[,9] <- rbinom(n=n_dev+n_val,1,.3) # comorbidity 1
    X[,10] <- rbinom(n=n_dev+n_val,1,.7) # comorbidity 2
    X[,11] <- rbinom(n=n_dev+n_val,1,.6) # comborbidity 3
    X[,12] <- rbinom(n=n_dev+n_val,1,.2) # comorbidity 4
    X[,13] <- rbinom(n_dev+n_val,1,.4) # comorbidity 5
    
    # generate 2 correlated count variables
    mu <- c(1,3)
    Sigma <- matrix(.7, nrow=2, ncol=2) + diag(2)*.3
    rawvars <- mvrnorm(n=n_dev+n_val,mu=mu, Sigma=Sigma)
    visits <- qpois(pnorm(rawvars),3)
    X[,14] <- visits[,1] # visits 1
    X[,15] <- visits[,2] # visits 2
    rm(mu, Sigma, rawvars, visits)
    
    # DGM that works for all versions
    xb1 <- -8.25 + .2*X[,1] + .24*((X[,7]*X[,10])) - .3*X[,3] + 
      .21*sqrt(X[,14]) - .9*X[,9] + .9*X[,11] + .01*sin(X[,5])
    
    xb2 <- -1.95 + .04*X[,1] + .5*((X[,7]*X[,10])) - .03*X[,3] + 
      .032*sqrt(X[,14]) - .02*X[,9] + .003*X[,11] + .31*sin(X[,5])
    xb1 <- xb1*1.8
    xb2 <- xb2*1.8
    p1 <- exp(xb1)/(1+exp(xb1)+exp(xb2))
    p2 <- exp(xb2)/(1+exp(xb1)+exp(xb2))
    p3 <- 1/(1+exp(xb1)+exp(xb2))
    pmat <- matrix(data=c(p1, p2,p3),nrow=n_dev+n_val,ncol=3)
    cmat <- colnames(as.data.frame(pmat))[max.col(as.data.frame(pmat),ties.method="first")] # convert probabilities to classes based on max pr
    Y <- case_when(cmat=="V1" ~ 1,
                   cmat=="V2" ~ 2,
                   TRUE ~ 3)
    xs <- c(8,9,10,12,2,4,5,14,15,6,13)
    
    dat <- as.data.frame(cbind(Y, X[,c(xs)]))
    colnames(dat) <- c("Y",paste0("X",xs))
    dat <- dat %>% mutate(id=row_number())
    # make survival data -- right censoring, Weibull dist.
    weib_shape <- 1
    weib_scale <- 90
    unif_weib <- runif(n_dev+n_val,0,1)
    Y_eff <- c(0,-1,1.5) # effect of class on survival
    xb_surv <- dat$X6 + Y_eff[dat$Y] 
    event <- ceiling((-log(unif_weib)*weib_scale*exp(-xb_surv))^(1/weib_shape))
    # censoring
    obsLen <- 365 #  1yr
    cens <- rep(obsLen, n_dev+n_val)
    dat$obsTime <- pmin(event, cens)
    dat$status <- event <=cens
    # use caret for stratified split into dev/val
    dat_split <- createFolds(factor(dat$Y), k=2, list=T) # not perfectly even b/c of stratification on y
    dat <- dat %>% mutate(foldid=ifelse((id %in% dat_split[[1]]),1,2))
    dat_dev <- dat %>% filter(foldid==1) %>% select(-foldid)
    dat_val <- dat %>% filter(foldid==2) %>% select(-foldid)
    rm(dat_split)
    
    # use caret for stratified splits into I1/I2 of dev data for weighted appraoch
    dat_split <- createFolds(factor(dat_dev$Y),k=2,list=T)
    # cant do by IDs here b/c have already split data, do by index instead
    split_index <- dat_split[[1]]
    dat_dev_I1 <- dat_dev[dat_split[[1]],]
    dat_dev_I2 <- dat_dev[dat_split[[2]],]
    
    # Naive: Use development data for algorithm fit, validation data for prediction ##
    c_drops <- c("id","X16", "obsTime","status") 
    dat_dev_fit <- dat_dev[ , !(names(dat_dev) %in% c_drops)] 
    # Weighted: use development data foldid==1 for alg fit
    dat_dev_I1_fit <- dat_dev_I1[, !names(dat_dev_I1) %in% c_drops]

    # just do one classification alg: multinomial logistic reg
    fit_naive <- nnet::multinom(Y ~ ., data=dat_dev_fit, maxit=500)
    fit_wtd <- nnet::multinom(Y ~ ., data=dat_dev_I1_fit, maxit=500)

    # naive sample (non-BS) prediction
    pred_val <- predict(fit_naive, dat_val, type="prob")
    y_pred_val <- max.col(pred_val) 
    pred_val <- as.data.frame(cbind(pred_val, y_pred_val, dat_val$id, dat_val$Y))
    colnames(pred_val) <- c("p1", "p2", "p3","y_pred", "id", "y_obs") 

    # weighted sample (non-bs) prediction in dev I2 -- for labeling
    pred_I2 <- predict(fit_wtd, dat_dev_I2, type="prob")
    y_max_I2 <- max.col(pred_I2)
    pred_I2 <- as.data.frame(cbind(pred_I2,y_max_I2, dat_dev_I2$id, dat_dev_I2$Y))
    colnames(pred_I2) <- c("p1", "p2", "p3","y_max", "id", "y_obs") 
    rm(y_max_I2)
    
    # naive coverage
    pred_tmp <- pred_val %>% transmute(g1=ifelse(y_pred==1,1,0),
                                g2=ifelse(y_pred==2,1,0),
                                g3=ifelse(y_pred==3,1,0))
    sample_coverage <- sapply(1:3,function(k)mean(pred_tmp[k==pred_val[,"y_obs"],k]))
    
    # naive classification performance
    sample_perf <- perfFun(pred_val)
    rm(pred_tmp)
    simrep_sample_perf_naive[[r]] <- sample_perf
    ## LABEL Function ##
    split_label_results <- splitLABELFun(data=pred_I2)

    ## Wtd Predictions in Validation Data ##
    pred_val_wtd <- predict(fit_wtd, dat_val, type="prob") 
    y_max_val_wtd <- max.col(pred_val_wtd[,1:3])
    pred_val_wtd <- cbind(pred_val_wtd, y_max_val_wtd, dat_val$Y)
    colnames(pred_val_wtd) <- c("p1", "p2", "p3", "y_pred","y_obs")
    rm(y_max_val_wtd)
    
    ## Label Predicted Classes in Validation Data based on thresholds set in Dev I2 ##
    label_results <- LABELFun(pred_val_wtd, split_label_results$thlds_est$class_thlds)
    print("label results checkpt")
    # no classification performance in single empirical sample for weighted BS method 
      # b/c no one single label assigned
    
    ## survival estimation -- observed Y
    km_val_obs <- survfit(Surv(obsTime, status) ~ Y, type="kaplan-meier",
                         data=dat_val, conf.type="plain")
    
    # observed survival probabilties at certain time points
    timeFun <- function(sample,time){
      surv <- summary(sample, times=time, extend=T)$surv
      lb <- summary(sample, times=time, extend=T)$lower
      ub <- summary(sample, times=time, extend=T)$upper
      out <- list(surv=surv,lb=lb, ub=ub)
      return(out)}
    times <- c(30,90,365) 
    sample_times_obs <- lapply(times, timeFun, sample=km_val_obs)
    names(sample_times_obs) <- paste0("d",times)
    
    # median survival time
    sample_50_obs <- quantile(km_val_obs, probs=.5, conf.int=T)
    simrep_km_obs[[r]] <- list(sample_times_obs=sample_times_obs, 
                               sample_median_obs=sample_50_obs)

    ## naive: survival estimation -- predicted Y
    pred_val <- dat_val %>% select(id, obsTime,status) %>% 
        left_join(pred_val,.,by="id") %>% 
      mutate(y_pred_fac=factor(y_pred,levels=c("1","2","3"), ordered=T)) 
    km_val_pred <- survfit(Surv(obsTime, status) ~ y_pred_fac,
                           type="kaplan-meier",
                           data=pred_val, conf.type="plain")
    sample_times_pred <- lapply(times, timeFun, sample=km_val_pred)
    names(sample_times_pred) <- paste0("d", times)
    sample_50_pred <- quantile(km_val_pred, probs=.5, conf.int=T)
    simrep_km_val_pred_sample[[r]] <- list(sample_times_pred=sample_times_pred,
                                           sample_median_pred=sample_50_pred)

    # data for bootstrap
    dat_bs <- cbind(as.data.frame(pred_val),label_results$Hstar_classcov)
    
    ## Bootstrap Function: To get variation around class coverage,
      # classification performance, and survival outcomes based on naive procedure ## 
    bsFun <- function(data){
      # naive 
      Hstar_naive <- data %>% transmute(g1=ifelse(y_pred==1,1,0),
                               g2=ifelse(y_pred==2,1,0),
                               g3=ifelse(y_pred==3,1,0))
      
      # coverage -- just to get BS-based CIs
      coverage_class_naive <- sapply(1:3,function(k)mean(Hstar_naive[k==data["y_obs"],k]))
      
      # km estimation
      # observed survival by observed stage for each resample
      km_obs <- survfit(Surv(obsTime,status)~ y_obs,type="kaplan-meier",
                        data=data,conf.type="plain")
      
      km_pred_naive <- survfit(Surv(obsTime, status) ~ y_pred_fac,type="kaplan-meier",
                                            data=data, conf.type="plain")

      out_naive <- list(coverage=coverage_class_naive, label_cols=data,
                  km_pred=km_pred_naive, obs=km_obs)
      
      # wtd
      Hstar_wtd <- data[,c("g1", "g2", "g3")]
      
      # coverage -- just to get BS-based CIs
      coverage_class_wtd <- sapply(1:3, function(k)mean(Hstar_wtd[k==data["y_obs"],k]))
      # identify ambiguous Hstars
      ambFun <- function(x){case_when(
        rowSums(x[,c("g1", "g2", "g3")])==0 ~ "0",
        rowSums(x[,c("g1", "g2", "g3")])==1 ~ "1",
        rowSums(x[,c("g1", "g2", "g3")])==2 ~ "2",
        TRUE ~ "3")}
      label_data <- data %>% mutate(amb_tmp=ambFun(.),
                                    amb_flag=factor(amb_tmp,
                                                    levels=c("0","1","2","3"), ordered=T)) %>% select(-amb_tmp)
      
      #######**  ASSIGNMENT OF LABEL **########
      # identify predicted class to use: if single label, seleted with pr==1
      # if 2 labels, select one fo the labels with pr==.5
      # if 0 or 3 labels, select one of labels with pr==.33 
      # could also assigned based on max pr for 0 (or any of them)
      labelFun <- function(x){case_when(
        x[,"amb_flag"]==0 ~ max.col(x[,c("g1", "g2", "g3")], ties.method="random"), # randomly select class 1-3 with equal probability
        # another option for 0 is to select class w/ max prob, but keeping simple for now
        x[,"amb_flag"]==3 ~ max.col(x[,c("g1", "g2", "g3")], ties.method="random"), # randomly select class 1-3 with equal probability
        x[,"amb_flag"]==2 ~ max.col(x[,c("g1", "g2", "g3")], ties.method="random"), # randomly select betwwen two assigned labels
        TRUE ~ max.col(x[,c("g1", "g2", "g3")]) # only one class assigned
      )}
      # max.col selects random col that is one of maxima
      label_data <- label_data %>% mutate(y_pred_wtd=labelFun(.))
      
      label_data <- label_data %>% mutate(y_pred_fac_wtd=factor(y_pred_wtd,
                                                            levels=c("1","2","3"), ordered=T))
      
      keep_cols <- label_data %>% select(amb_flag, g1, g2, g3,y_pred, y_pred_fac,
                                         y_pred_wtd, y_obs, p1, p2, p3, y_pred_fac_wtd)
      
      km_pred_wtd <- survfit(Surv(obsTime, status) ~ y_pred_fac_wtd,type="kaplan-meier",
                         data=label_data, conf.type="plain")
      
      out_wtd <- list(coverage=coverage_class_wtd, label_cols=keep_cols,
                  km_pred=km_pred_wtd)
      
      out <- list(naive=out_naive, wtd=out_wtd)
      return(out)
    }
    
    ## Generate Resamples for Bootstrap ## 
    resamples <- lapply(1:n_bs, function(i)dat_bs[sample(nrow(dat_bs),n_val,replace=T),])
    print("resamples checkpt")
    ## Apply BS function to Resamples ##
    bs_out <- lapply(resamples, bsFun)
    print("bs out checkpt")

    ## Measures that vary across BS iterations ##
    # coverage
    simrep_class_coverage_naive[[r]] <- round(covFun(lapply(lapply(bs_out, `[[`, "naive"), `[[`, "coverage"), sample_coverage), 3)[1,]
    simrep_class_coverage_wtd[[r]] <- round(covFun(lapply(lapply(bs_out, `[[`, "wtd"), `[[`, "coverage"), sample_cov=label_results$coverage_class), 3)
    print("cov checkpt")
    
    # ambiguity
    simrep_ambiguity[[r]] <- ambFun(lapply(bs_out, `[[`, "wtd"))
    
    # classification performance 
    simrep_class_measures_naive_pct[[r]]  <- perfBSFun_pct(lapply(bs_out, `[[`, "naive"))
    simrep_class_measures_wtd[[r]] <- perfWTDBSFun(lapply(bs_out, `[[`, "wtd"))
  
    # survival estimation + bias
    data_obs_bias <- lapply(lapply(bs_out,`[[`,"naive"),`[[`,"obs")
    simrep_km_naive_pct[[r]] <- medWTDBSFun(lapply(lapply(bs_out, `[[`,"naive"),`[[`, "km_pred"),data_obs_bias)
    simrep_km_wtd[[r]] <- medWTDBSFun(lapply(lapply(bs_out, `[[`,"wtd"),`[[`, "km_pred"),data_obs_bias)
    
    rm(bs_out, bsFun, c_drops, cens, cmat, dat, dat_bs,
       dat_dev, dat_dev_fit,dat_dev_I1, dat_dev_I1_fit, dat_dev_I2,
       dat_split, dat_val, event, fit_naive, fit_wtd,
       km_val_obs, km_val_pred, label_results, obsLen, p1, p2, p3,
       pmat,pred_I2, pred_val,pred_val_wtd, resamples, sample_50_obs,
       sample_coverage, sample_perf, sample_times_obs, split_label_results,
       times, unif_weib, weib_scale, weib_shape, X, xb_surv, xb1, 
       xb2,xs, Y, Y_eff, y_pred_val,data_obs_bias)
    
    print(paste0(r, " simrep done"))
    }
)


print("done with sims")


savepath <- c("~/sim_output/")

saveRDS(simrep_km_naive_pct, paste0(savepath, "simrep_km_naive_pct_", repname, ".RDS"))
saveRDS(simrep_km_wtd, paste0(savepath, "simrep_km_wtd_", repname, ".RDS"))
saveRDS(simrep_km_obs, paste0(savepath, "simrep_km_obs_", repname, ".RDS"))
saveRDS(simrep_km_val_pred_sample, paste0(savepath, "simrep_km_val_pred_sample_", repname, ".RDS"))

saveRDS(simrep_thresholds, paste0(savepath, "simrep_thresholds_", repname, ".RDS"))
saveRDS(simrep_ambiguity, paste0(savepath, "simrep_ambiguity_", repname, ".RDS"))
saveRDS(simrep_class_coverage_naive, paste0(savepath, "simrep_class_coverage_naive_", repname, ".RDS"))
saveRDS(simrep_class_coverage_wtd, paste0(savepath, "simrep_class_coverage_wtd_", repname, ".RDS"))

saveRDS(simrep_class_measures_naive_pct, paste0(savepath, "simrep_class_measures_naive_pct_", repname, ".RDS"))
saveRDS(simrep_class_measures_wtd, paste0(savepath, "simrep_class_measures_wtd_", repname, ".RDS"))

saveRDS(simrep_sample_perf_naive, paste0(savepath, "simrep_sample_perf_naive_", repname, ".RDS"))




