#########################################################
# CancerCLAS Multiclass Prediction
# Empirical Data
# Naive No Variation, Naive Bootstrap, Weighted Bootstrap
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


# source misc_funs.R: contains function for calculating class measures
source("~/misc_funs_empirical_v3.R")

# Generate some example training data
set.seed(33)

#### Data Set Up ####
# select number of bootstrap iterations
n_bs <- 500
# alpha for setting label thresholds
alpha <- .10

### Development data: 2010-2011 ###
load("~/lung_20102011_072819.RData")
data_dev <- lung_20102011_072819; rm(lung_20102011_072819)
order_dev <- colnames(data_dev) # column order -- apply to validation data below
# shuffle data
data_dev <- data_dev[sample(nrow(data_dev)),]
# assign row ids for splitting
data_dev <- data_dev %>% mutate(id=row_number())
# use caret for stratified CV folds
data_split <- createFolds(factor(data_dev$StageGroup_AJCC6), k=2, list=T)
# not perfectly even b/c of stratification on y
# then make foldids a column
data_dev <- data_dev %>% mutate(foldid = ifelse((id %in% data_split[[1]]),1,2))
foldid <- data_dev$foldids
# remove ref groups, survival, and others for fitting
data_dev_fit <- data_dev %>% dplyr::select(-Early_Stage, -StageGroup_Predicted, 
                               -Region_West, -Surgery_Segmental, 
                               -CauseOfDeath_LungCancer, 
                               -Days_First_Chemo_to_Death, -Censored)
rm(data_split)
id_drops <- c("id", "foldid")

## Naive fit on development data ##
fitMN <- function(data){
  data <- data[, !(names(data) %in% id_drops)]
  fit_mn <- nnet::multinom(StageGroup_AJCC6 ~., data=data, maxit=500)
  return(fit_mn)}
fitGLMNET <- function(data,alpha){
  x_drops <- c("id", "foldid", "StageGroup_AJCC6")
  xdata <- as.matrix(data[, !names(data) %in% x_drops])
  fit_glmnet <- glmnet(x=xdata, y=as.factor(data[,"StageGroup_AJCC6"]),family="multinomial", alpha=alpha)
  return(fit_glmnet)}
fitRF <- function(data){
  data <- data[, !(names(data) %in% id_drops)]
  fit_rf <- randomForest(as.factor(StageGroup_AJCC6) ~., data=data, 
                         ntree=500, nodesize=250,
                         strata=as.factor(data$StageGroup_AJCC6))
  return(fit_rf)}
fitGAM <- function(data){
  # make Y be numeric, scaled from 0-2
  data[,"StageGroup_AJCC6"] <- as.numeric(data$StageGroup_AJCC6)-1
  # create gam formula (from misc_funs); cubic splines with k=3 knots
  # could select number of knots via internal cross validation, but 
  # because we have so many variables this would take a long time, 
  # so going for practical approach of setting uniform number
  # same thing with the smoothing penalty
  f <- CreateGAMFormula(data=data[,!names(data) %in% id_drops],
                        y="StageGroup_AJCC6", type="regspline")
  # remove id and foldid in formula step
  
  f1 <- f[[1]]
  f2 <- f[[2]]
  # s=0.6 and k=3 for all terms
  fit <- mgcv::gam(list(f1,f2), data=data,family=multinom(K=2))
  fit_gam <- list(fit=fit, data=data)
  return(fit_gam)
}
fitXGB <- function(data){
  data <- data[,!names(data) %in% id_drops]
  data[,"StageGroup_AJCC6"] <- as.numeric(data$StageGroup_AJCC6)-1
  data <- as.matrix(data) 
  
  fit <- xgboost(data=subset(data, select=-StageGroup_AJCC6),
                 label=data[,"StageGroup_AJCC6"],
                 max.depth=3, eta=1, nthread=3, nrounds=2,
                 objective="multi:softprob", num_class=3,
                 eval_metric="mlogloss")
  fit_xgb <- list(fit=fit, data=data)
  return(fit_xgb)
}

devNaiveFitFun <- function(data){
  # mn
  fit_mn <- fitMN(data=data)
  # lasso
  fit_lasso <- fitGLMNET(data=data,alpha=1)
  # ridge
  fit_ridge <- fitGLMNET(data=data, alpha=0)
  # enet 
  fit_enet <- fitGLMNET(data=data, alpha=.5)
  # rf
  fit_rf  <- fitRF(data)
  # gam
  fit_gam <- fitGAM(data)
  # xgb
  fit_xgb <- fitXGB(data)
  
  out <- list(fit_mn=fit_mn, fit_lasso=fit_lasso, fit_ridge=fit_ridge,
              fit_enet=fit_enet, fit_rf=fit_rf, fit_gam=fit_gam, fit_xgb=fit_xgb)
}
dev_naive_results <- devNaiveFitFun(data=data_dev_fit)


### now use split for algorithm fitting
devWtdFitFun <- function(data,thld_fold_num){
  holdoutIndex <- which(data[,"foldid"]==thld_fold_num, arr.ind=T)
  holdoutFold <- data[holdoutIndex,]
  fitFolds <- data[-holdoutIndex,]
  
  # multinomial algorithm   
  holdoutFold_mn <- as.matrix(holdoutFold[, !names(holdoutFold) %in% id_drops])
  fit_mn <- fitMN(data=fitFolds)
  pred_mn <- predict(fit_mn, holdoutFold_mn, type="prob")
  y_pred <- max.col(pred_mn)
  pred_mn <- cbind(pred_mn, y_pred, holdoutFold$id, holdoutFold$foldid, holdoutFold$StageGroup_AJCC6)
  colnames(pred_mn) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  
  # lasso -- specify variable thresholds
  #drop vars and convert to matrix for input
  x_drops <- c("id", "foldid", "StageGroup_AJCC6")
  fitFolds_glmnet <- as.matrix(fitFolds[, !names(fitFolds) %in% x_drops])
  # use fitFolds_glmnet for prediction below
  holdoutFold_glmnet <- as.matrix(holdoutFold[, !names(holdoutFold) %in% x_drops])
  fit_lasso <- fitGLMNET(data=data,alpha=1) 
  
  # get the minimum lambda, get coefs for each category
  lambda_min <- min(fit_lasso$lambda)
  fit_lasso_coefs <- coef(fit_lasso, s = lambda_min)
  pred_lasso <- predict(fit_lasso, as.matrix(holdoutFold_glmnet),
                            type="response", s=lambda_min)
  fit_lasso_values <- predict(fit_lasso, as.matrix(fitFolds_glmnet),
                              type="response", s=lambda_min)
  y_pred <- max.col(pred_lasso[,1:3,])
  pred_lasso <- as.matrix(cbind(as.data.frame(pred_lasso),y_pred,
                                    holdoutFold$id, holdoutFold$foldid, 
                                    holdoutFold$StageGroup_AJCC6))
  colnames(pred_lasso) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  print("lasso")
  
  # ridge
  fit_ridge <- fitGLMNET(data=data,alpha=0)
  lambda_min_ridge <- min(fit_ridge$lambda)
  pred_ridge <- predict(fit_ridge, as.matrix(holdoutFold_glmnet),
                        type="response", s=lambda_min_ridge)
  fit_ridge_values <- predict(fit_ridge, as.matrix(fitFolds_glmnet),
                              type="response", s=lambda_min_ridge)
  y_pred <- max.col(pred_ridge[,1:3,])
  pred_ridge <- as.matrix(cbind(as.data.frame(pred_ridge),y_pred,
                                holdoutFold$id, holdoutFold$foldid, 
                                holdoutFold$StageGroup_AJCC6))
  colnames(pred_ridge) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  print("ridge")
  
  # elastic net
  fit_enet <- fitGLMNET(data=data, alpha=.5)
  lambda_min_enet <- min(fit_enet$lambda)
  pred_enet <- predict(fit_enet, as.matrix(holdoutFold_glmnet),
                       type="response", s=lambda_min_enet)
  fit_enet_values <- predict(fit_enet, as.matrix(fitFolds_glmnet),
                             type="response", s=lambda_min_enet)
  y_pred <- max.col(pred_enet[,1:3,])
  pred_enet <- as.matrix(cbind(as.data.frame(pred_enet), y_pred,
                               holdoutFold$id, holdoutFold$foldid, 
                               holdoutFold$StageGroup_AJCC6))
  colnames(pred_enet) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  print("enet")
  # random forest
  holdoutFold_rf <- holdoutFold[,!(names(holdoutFold) %in% id_drops)]
  fit_rf <- fitRF(fitFolds)
  pred_rf <- predict(fit_rf, holdoutFold_rf, type="prob")
  y_pred <- max.col(pred_rf)
  pred_rf <- cbind(pred_rf,y_pred, holdoutFold$id, holdoutFold$foldid, 
                   holdoutFold$StageGroup_AJCC6)
  colnames(pred_rf) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  print("rf")
  
  # GAM
  fit_gam <- fitGAM(data=fitFolds)
  # prediction step
  holdoutFold_gam <- holdoutFold
  holdoutFold_gam[,"StageGroup_AJCC6"] <- as.numeric(holdoutFold_gam[,"StageGroup_AJCC6"])-1
  pred_gam <- predict(fit_gam$fit, newdata=holdoutFold_gam, type="response")
  fit_gam_values <- predict(fit_gam$fit, newdata=fit_gam$data, type="response")
  y_pred <- max.col(pred_gam)
  pred_gam <- cbind(pred_gam,y_pred, holdoutFold$id, holdoutFold$foldid, 
                    holdoutFold$StageGroup_AJCC6)
  colnames(pred_gam) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  print("gam")
  
  # Boosting
  # convert fitFolds and holdoutFold to matrices
  holdoutFold_xgb <- holdoutFold
  holdoutFold_xgb$StageGroup_AJCC6 <-  as.numeric(holdoutFold_xgb$StageGroup_AJCC6)-1
  holdoutFold_xgb_noids <- holdoutFold_xgb[,!names(holdoutFold_xgb) %in% id_drops]
  holdoutFold_xgb_noids <- as.matrix(holdoutFold_xgb_noids)
  holdoutFold_xgb_noids <- xgb.DMatrix(data=holdoutFold_xgb_noids[,-95])
  # remove label column (StageGroup_AJCC6)
   fit_xgb <- fitXGB(fitFolds)
  # predict outputs data*nclass vector, turn into
  # data*nclass matrix
  pred_xgb <- matrix(predict(fit_xgb$fit, holdoutFold_xgb_noids),
                     nrow=nrow(holdoutFold_xgb_noids), byrow=T)
  # havee to drop out StageGroup_AJCC6 from fitFolds_xgb before predicting
  fit_xgb_values <- matrix(predict(fit_xgb$fit, fit_xgb$data[,-95]),
                           nrow=nrow(fit_xgb$data), byrow=T)
  y_pred <- max.col(pred_xgb)
  pred_xgb <- cbind(pred_xgb,y_pred, holdoutFold$id, holdoutFold$foldid, 
                    holdoutFold$StageGroup_AJCC6)
  colnames(pred_xgb) <- c("g1", "g2", "g3","y_pred", "id", "foldid", "y_obs")
  
  
  out <- list(fit_mn=fit_mn, pred_mn=pred_mn,
              fit_lasso = fit_lasso, pred_lasso=pred_lasso,
              fit_lasso_coefs= fit_lasso_coefs,fit_lasso_values=fit_lasso_values,
              fit_ridge=fit_ridge, pred_ridge=pred_ridge,fit_ridge_values=fit_ridge_values,
              fit_enet=fit_enet, pred_enet=pred_enet,fit_enet_values=fit_enet_values,
              fit_rf=fit_rf, pred_rf=pred_rf,
              fit_gam=fit_gam, pred_gam=pred_gam,fit_gam_values=fit_gam_values,
              fit_xgb=fit_xgb, pred_xgb=pred_xgb, fit_xgb_values=fit_xgb_values)
  return(out)
}
# fit on I1 fold, predict on I2
dev_wtd_results <- devWtdFitFun(data=data_dev_fit, thld_fold_num=2)

# define thresholds and create label sets for I2
devLABELFun <- function(data,dev_results,...){
  
  # threshold estimation
  conf_thlds <- function(phat,Y2,class_alphas){
    m = nrow(phat)
    K = 3
    score = rep(0,m)
    for(i in 1:m){
      score[i] = phat[i,Y2[i]]
    }

    ## Class-specific coverage
    class_thlds <- rep(NA,K)
    for (k in 1:3){
      class_thlds[k] <- sort(score[Y2==k])[ceiling(class_alphas[k]*(sum(Y2==k)+1)-1)]
    }
    
    return(class_thlds=class_thlds)
    
  } 
  
  # use fitted values and obs Ys to estimate thresholds, and then apply to predicted values
  thlds_est_mn <- conf_thlds(phat=dev_results$pred_mn, 
                             Y2=data[,"StageGroup_AJCC6"],class_alphas=c(.1,.1,.1))
  # do this for all the other algs
  thlds_est_lasso <- conf_thlds(phat=as.data.frame(dev_results$pred_lasso), 
                                Y2=data[,"StageGroup_AJCC6"], class_alphas=c(.1,.1,.1))
  
  thlds_est_ridge <- conf_thlds(phat=as.data.frame(dev_results$pred_ridge),
                                Y2=data[,"StageGroup_AJCC6"], class_alphas=c(.1,.1,.1))
  
  thlds_est_enet <- conf_thlds(phat=as.data.frame(dev_results$pred_enet),
                               Y2=data[,"StageGroup_AJCC6"], class_alphas=c(.1,.1,.1))
  
  thlds_est_rf <- conf_thlds(phat=dev_results$pred_rf, 
                             Y2=data[,"StageGroup_AJCC6"],class_alphas=c(.1,.1,.1))
  
  thlds_est_gam <- conf_thlds(phat=dev_results$pred_gam, 
                              Y2=data[,"StageGroup_AJCC6"], class_alphas=c(.1,.1,.1))
  
  thlds_est_xgb <- conf_thlds(phat=dev_results$pred_xgb, 
                              Y2=data[,"StageGroup_AJCC6"], class_alphas=c(.1,.1,.1))
  
  hstarClassFun <- function(pred_name, thlds_name){
    tmp <- dev_results[[pred_name]]
    tmp <- tmp[,1:3]
    t(apply(tmp, 1,
            function(x){as.numeric(x >= thlds_name)}))
  }
  
  
  Hstar_classcov_mn <- hstarClassFun("pred_mn", thlds_est_mn)
  Hstar_classcov_lasso <- hstarClassFun("pred_lasso", thlds_est_lasso)
  Hstar_classcov_ridge <- hstarClassFun("pred_ridge", thlds_est_ridge)
  Hstar_classcov_enet <- hstarClassFun("pred_enet", thlds_est_enet)
  Hstar_classcov_rf <- hstarClassFun("pred_rf", thlds_est_rf)
  Hstar_classcov_gam <- hstarClassFun("pred_gam", thlds_est_gam)
  Hstar_classcov_xgb <- hstarClassFun("pred_xgb", thlds_est_xgb)
  
  out <- list(thlds_est_mn=thlds_est_mn,
              thlds_est_lasso=thlds_est_lasso, 
              thlds_est_ridge=thlds_est_ridge, 
              thlds_est_enet=thlds_est_enet, 
              thlds_est_rf=thlds_est_rf,
              thlds_est_gam=thlds_est_gam,  
              thlds_est_xgb=thlds_est_xgb)
  return(out)
}
data_dev_I2 <- data_dev %>% filter(foldid==2) %>% as.data.frame
data_dev_thlds <- devLABELFun(data=data_dev_I2, dev_results=dev_wtd_results)

### Validation data: 2012-2013 ###
load("~/lung_20122013_072819.RData")
data_val <- lung_20122013_072819; rm(lung_20122013_072819)
# make column names match 2010-2011 data/ arrange in same order
data_val <- data_val %>% select(order_dev)
# drop if Days_First_Chemo_to_Death is <0
data_val <- data_val %>% filter(Days_First_Chemo_to_Death>=0|is.na(Days_First_Chemo_to_Death))
# assign ids to keep track across resamples
data_val <- data_val %>% mutate(id=row_number())
data_val_pred <- data_val %>% select(-StageGroup_Predicted, -Early_Stage, -Region_West,
                        -Surgery_Segmental)


## WTD Validation Prediction ##
# use conditional probability estimators fit on I1_dev; thresholds based on I2_dev
# pull out fitted algorithms (fit on I1)
valWTDPredFun <- function(data, data_for_pred){
  # remove non-prediction vars from prediction for glmnet-based algs
  x_drops <- c("StageGroup_AJCC6", "CauseOfDeath_LungCancer", "Censored", "Days_First_Chemo_to_Death",
               "id", "foldid")
  data_glmnet <- as.matrix(data_val_pred[, !names(data_val_pred) %in% x_drops]) 
  
  # list to hold predication output
  out_list <- vector("list", 7) 
  
  ## predictions in 2012-2013 validation data ##
  # create matrix with predicted probabilities for each class, the 
  # single predicted class based on highest probability, and
  # observed class
  
  # mn
  pred_mn <- predict(dev_wtd_results$fit_mn, data_for_pred, type="prob") 
  y_pred <- max.col(pred_mn[,1:3])
  pred_mn <- cbind(pred_mn, y_pred, data$StageGroup_AJCC6, data$id)
  colnames(pred_mn) <- c("p1", "p2", "p3", "y_pred","y_obs", "id")
  out_list[[1]] <- as.data.frame(pred_mn)
  
  ## lasso 
  # pull out lambda min first (for bootstrap version, do this outside bootstrap)
  lambda_min <- min(dev_wtd_results$fit_lasso$lambda)
  pred_lasso <- predict(dev_wtd_results$fit_lasso, as.matrix(data_glmnet),
                        type="response", s=lambda_min)
  y_pred <- max.col(pred_lasso[,1:3,])
  pred_lasso <- as.matrix(cbind(as.data.frame(pred_lasso),y_pred, data$StageGroup_AJCC6, data$id))
  colnames(pred_lasso) <- c("p1", "p2", "p3","y_pred", "y_obs", "id")
  out_list[[2]] <- as.data.frame(pred_lasso)
  
  ## ridge
  lambda_min_ridge <- min(dev_wtd_results$fit_ridge$lambda)
  pred_ridge <- predict(dev_wtd_results$fit_ridge, as.matrix(data_glmnet),
                        type="response", s=lambda_min_ridge)
  y_pred <- max.col(pred_ridge[,1:3,])
  pred_ridge <- as.matrix(cbind(as.data.frame(pred_ridge), y_pred,data$StageGroup_AJCC6, data$id))
  colnames(pred_ridge) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[3]] <- as.data.frame(pred_ridge)
  
  ## enet
  lambda_min_enet <- min(dev_wtd_results$fit_enet$lambda)
  pred_enet <- predict(dev_wtd_results$fit_enet, as.matrix(data_glmnet),
                       type="response", s=lambda_min_enet)
  y_pred <- max.col(pred_enet[,1:3,])
  pred_enet <- as.matrix(cbind(as.data.frame(pred_enet), y_pred,data$StageGroup_AJCC6, data$id))
  colnames(pred_enet) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[4]] <- as.data.frame(pred_enet)
  
  ## rf
  pred_rf <- predict(dev_wtd_results$fit_rf, data_for_pred, type="prob")
  y_pred <- max.col(pred_rf[,1:3])
  pred_rf <- cbind(pred_rf,y_pred, data$StageGroup_AJCC6, data$id)
  colnames(pred_rf) <- c("p1", "p2", "p3","y_pred","y_obs","id")
  out_list[[5]] <- as.data.frame(pred_rf)
  
  ## gam
  pred_gam <- predict(dev_wtd_results$fit_gam$fit, newdata=data_for_pred, type="response")
  y_pred <- max.col(pred_gam[,1:3])
  pred_gam <- cbind(pred_gam,y_pred, data$StageGroup_AJCC6, data$id)
  colnames(pred_gam) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[6]] <- as.data.frame(pred_gam)
  
  ## xgb
  data_xgb <- as.matrix(data_for_pred)
  data_xgb_pred <- xgb.DMatrix(data=subset(data_xgb, select=c(-StageGroup_AJCC6, 
                                                              -CauseOfDeath_LungCancer,
                                                              -Censored,
                                                              -Days_First_Chemo_to_Death,
                                                              -id)),
                               label=data_xgb[,"StageGroup_AJCC6"])
  pred_xgb <- matrix(predict(dev_wtd_results$fit_xgb$fit, data_xgb_pred),
                     nrow=nrow(data_xgb_pred), byrow=T)
  y_pred <- max.col(pred_xgb[,1:3])
  pred_xgb <- cbind(pred_xgb,y_pred, data$StageGroup_AJCC6, data$id)
  colnames(pred_xgb) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[7]] <- as.data.frame(pred_xgb)
  
  return(out_list)
}
val_wtd_results <- valWTDPredFun(data=data_val, data_for_pred=data_val_pred)
alg_names <- c("mn", "lasso", "ridge", "enet", "rf", "gam","xgb")
names(val_wtd_results) <- alg_names
# function labeling predicted classes in 2012-2013 data based on thresholds set in 
# split 2010-2011 data
val_label_results <- LABELFun(val_wtd_results, thlds=data_dev_thlds, data_val=data_val)
# val label results captures wtd sample coverage, 
## Wtd Validation ambiguity ##
# val_label_results already captures wtd sample coverage, ambiguity

## Naive Validation Prediction ##
valNaivePredFun <- function(data, data_for_pred){
  # remove non-prediction vars from prediction for glmnet-based algs
  x_drops <- c("StageGroup_AJCC6", "CauseOfDeath_LungCancer", "Censored", "Days_First_Chemo_to_Death",
               "id", "foldid")
  data_glmnet <- as.matrix(data_val_pred[, !names(data_val_pred) %in% x_drops]) 
  
  # list to hold predication output
  out_list <- vector("list", 7) 
  
  ## predictions in 2012-2013 validation data ##
  # create matrix with predicted probabilities for each class, the 
  # single predicted class based on highest probability, and
  # observed class, and id
  
  # mn
  pred_mn <- predict(dev_naive_results$fit_mn, data_for_pred, type="prob") 
  y_pred <- max.col(pred_mn[,1:3])
  pred_mn <- cbind(pred_mn, y_pred, data$StageGroup_AJCC6,data$id)
  colnames(pred_mn) <- c("p1", "p2", "p3", "y_pred","y_obs","id")
  out_list[[1]] <- as.data.frame(pred_mn)
  
  ## lasso 
  # pull out lambda min first (for bootstrap version, do this outside bootstrap)
  lambda_min <- min(dev_naive_results$fit_lasso$lambda)
  pred_lasso <- predict(dev_naive_results$fit_lasso, as.matrix(data_glmnet),
                        type="response", s=lambda_min)
  y_pred <- max.col(pred_lasso[,1:3,])
  pred_lasso <- as.matrix(cbind(as.data.frame(pred_lasso),y_pred, data$StageGroup_AJCC6,data$id))
  colnames(pred_lasso) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[2]] <- as.data.frame(pred_lasso)
  
  ## ridge
  lambda_min_ridge <- min(dev_naive_results$fit_ridge$lambda)
  pred_ridge <- predict(dev_naive_results$fit_ridge, as.matrix(data_glmnet),
                        type="response", s=lambda_min_ridge)
  y_pred <- max.col(pred_ridge[,1:3,])
  pred_ridge <- as.matrix(cbind(as.data.frame(pred_ridge), y_pred,data$StageGroup_AJCC6,data$id))
  colnames(pred_ridge) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[3]] <- as.data.frame(pred_ridge)
  
  ## enet
  lambda_min_enet <- min(dev_naive_results$fit_enet$lambda)
  pred_enet <- predict(dev_naive_results$fit_enet, as.matrix(data_glmnet),
                       type="response", s=lambda_min_enet)
  y_pred <- max.col(pred_enet[,1:3,])
  pred_enet <- as.matrix(cbind(as.data.frame(pred_enet), y_pred,data$StageGroup_AJCC6,data$id))
  colnames(pred_enet) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[4]] <- as.data.frame(pred_enet)
  
  ## rf
  pred_rf <- predict(dev_naive_results$fit_rf, data_for_pred, type="prob")
  y_pred <- max.col(pred_rf[,1:3])
  pred_rf <- cbind(pred_rf,y_pred, data$StageGroup_AJCC6,data$id)
  colnames(pred_rf) <- c("p1", "p2", "p3","y_pred","y_obs","id")
  out_list[[5]] <- as.data.frame(pred_rf)
  
  ## gam
  pred_gam <- predict(dev_naive_results$fit_gam$fit, newdata=data_for_pred, type="response")
  y_pred <- max.col(pred_gam[,1:3])
  pred_gam <- cbind(pred_gam,y_pred, data$StageGroup_AJCC6,data$id)
  colnames(pred_gam) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[6]] <- as.data.frame(pred_gam)
  
  ## xgb
  data_xgb <- as.matrix(data_for_pred)
  data_xgb_pred <- xgb.DMatrix(data=subset(data_xgb, select=c(-StageGroup_AJCC6, 
                                                              -CauseOfDeath_LungCancer,
                                                              -Censored,
                                                              -Days_First_Chemo_to_Death,
                                                              -id)),
                               label=data_xgb[,"StageGroup_AJCC6"])
  pred_xgb <- matrix(predict(dev_naive_results$fit_xgb$fit, data_xgb_pred),
                     nrow=nrow(data_xgb_pred), byrow=T)
  y_pred <- max.col(pred_xgb[,1:3])
  pred_xgb <- cbind(pred_xgb,y_pred, data$StageGroup_AJCC6,data$id)
  colnames(pred_xgb) <- c("p1", "p2", "p3","y_pred", "y_obs","id")
  out_list[[7]] <- as.data.frame(pred_xgb)
  
  return(out_list)
}
val_naive_results <- valNaivePredFun(data=data_val, data_for_pred=data_val_pred)
names(val_naive_results) <- alg_names

## Naive Validation Sample Coverage ##
covNaiveSampleFun <- function(data){
  pred_tmp <- data %>% transmute(g1=ifelse(y_pred==1,1,0),
                                              g2=ifelse(y_pred==2,1,0),
                                              g3=ifelse(y_pred==3,1,0))
  sample_coverage <- sapply(1:3,function(k)mean(pred_tmp[k==data[,"y_obs"],k]))
  return(sample_coverage)
}
val_naive_sample_coverage <- lapply(val_naive_results, covNaiveSampleFun)

## In sample naive no-boot classification performance
val_naive_pred_in_sample <- lapply(val_naive_results, perfFun)

## In sample km based on observed class -- same across all algsm just use data_val
# create censoring indicator (1=dead/not censored, 0=alive  to play well with Surv function)
eventFun <- function(x){ifelse(is.na(x[,"Days_First_Chemo_to_Death"]),0,1)}
timeFun <- function(x){ifelse(x[,"event"]==0,max(x[,"Days_First_Chemo_to_Death"],na.rm=T),
                              x[,"Days_First_Chemo_to_Death"])}
time2Fun <- function(x){ifelse(x[,"time"]==0,1,x[,"time"])}

data_val_tmp <- data_val %>% mutate(event=eventFun(.)) %>% mutate(time=timeFun(.)) %>%
  mutate(time2=time2Fun(.))

val_km_obs_in_sample <-  survfit(Surv(time2,event)~StageGroup_AJCC6,type="kaplan-meier",
                                 data=data_val_tmp,conf.type="plain")
# observed survival probabilties at certain time points
survTimeFun <- function(sample,time){
  surv <- summary(sample, times=time, extend=T)$surv
  lb <- summary(sample, times=time, extend=T)$lower
  ub <- summary(sample, times=time, extend=T)$upper
  out <- list(surv=surv,lb=lb, ub=ub)
  return(out)}
times <- c(90,365) 
sample_times_obs <- lapply(times, survTimeFun, sample=val_km_obs_in_sample)
names(sample_times_obs) <- paste0("d",times)
# median survival time -- ignore for stage I/II
sample_50_obs <- quantile(val_km_obs_in_sample, probs=.5, conf.int=T)
sample_km_obs <- list(sample_times_obs=sample_times_obs, sample_median_obs=sample_50_obs)
rm(sample_times_obs, sample_50_obs)
## In sample km based on predicted class (naive no-boot)
naiveKMPredFun <- function(data){
  tmp <- data_val_tmp %>% select(id, time2, event) %>% left_join(data,.,by="id") %>%
    mutate(y_pred_fac=factor(y_pred,levels=c("1","2","3"),ordered=T))
  km <- survfit(Surv(time2,event)~ y_pred_fac, type="kaplan-meier", data=tmp,
                conf.type="plain")
  sample_times <- lapply(times, survTimeFun, sample=km)
  names(sample_times) <- paste0("d",times)
  sample_50 <- quantile(km, probs=.5, conf.int=T, na.rm=T)
  out <- list(sample_times_pred=sample_times, sample_median_pred=sample_50)
  return(out)
}
sample_km_pred_naive <- lapply(val_naive_results, naiveKMPredFun)

### Bootstrap: Naive and Weighted ###
# keep only Hstar data.frames
Hstar_list <- val_label_results$Hstar
# join pred values, data, and labels together
  # note that joining in rest of data is unncessary if just doing KM
val_pred_labels_data <- map2(val_wtd_results, Hstar_list, ~cbind(.x, .y))
val_pred_labels_data <- map(val_pred_labels_data, ~left_join(.x,data_val,by=c("id")))


bsFun <- function(data){
  # create variables for survival estimation
  # create censoring indicator (1=dead/not censored, 0=alive  to play well with Surv function)

  data <- data %>% mutate(event=eventFun(.))
  data <- data %>% mutate(time=timeFun(.))
  data <- data %>% mutate(time2=time2Fun(.))
  
  # naive 
  Hstar_naive <- data %>% transmute(g1=ifelse(y_pred==1,1,0),
                                    g2=ifelse(y_pred==2,1,0),
                                    g3=ifelse(y_pred==3,1,0))
    
  coverage_class_naive <- sapply(1:3,function(k)mean(Hstar_naive[k==data["y_obs"],k]))
  

  # km estimation
  # observed survival by observed stage for each resample
  km_obs <- survfit(Surv(time2,event)~StageGroup_AJCC6,type="kaplan-meier",
                        data=data,conf.type="plain")
  # observed survival by predicted stage for each resample
  km_pred_naive <- survfit(Surv(time2,event)~y_pred,type="kaplan-meier",
                           data=data,conf.type="plain")
  out_naive <- list(coverage=coverage_class_naive,km_pred=km_pred_naive,obs=km_obs)
  
  # wtd 
  Hstar_wtd <- data[,c("g1", "g2", "g3")]
  
  # coverage
  coverage_class_wtd <- sapply(1:3,
                           function(k)mean(Hstar_wtd[k==data["y_obs"],k]))
  
  # identify ambiguous Hstars
  ambFun <- function(x){case_when(
    rowSums(x[,c("g1","g2","g3")])==0 ~"0",
    rowSums(x[,c("g1","g2","g3")])==1 ~"1",
    rowSums(x[,c("g1","g2","g3")])==2 ~"2",
    TRUE~"3")}
  label_data <- data %>% mutate(amb_tmp=ambFun(.),
                                amb_flag=factor(amb_tmp,levels=c("0","1","2","3"),ordered=T)) %>% select(-amb_tmp)
  
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
  keep_cols <- label_data %>% select(amb_flag, g1, g2, g3,y_pred,
                                     y_pred_wtd, y_obs, p1, p2, p3)
  
  km_pred_wtd <- survfit(Surv(time2,event)~y_pred_wtd,type="kaplan-meier",
                        data=label_data,conf.type="plain")
  
  out_wtd <- list(coverage=coverage_class_wtd,label_cols=keep_cols,
                  km_pred=km_pred_wtd)
  
  out <- list(naive=out_naive,wtd=out_wtd)
  return(out)
}

# making the resamples takes <1 min each; total resmaples list is 28 Gb
resamplesFun <- function(alg){
  tmp <- lapply(1:500, function(i)val_pred_labels_data[[alg]][sample(nrow(val_pred_labels_data[[alg]]),nrow(data_val),replace=T),])
  return(tmp)
}
resamples_list <- lapply(alg_names, resamplesFun)
names(resamples_list) <- alg_names

# apply bsFun across all alg resamples
system.time(bs_out_list <- lapply(resamples_list, function(x){lapply(x, bsFun)}))
rm(resamples_list)
## Measures that vary across BS iterations ##
bsMeasFun <- function(data, alg){
  # coverage - to get bs-based CIs
  cov_naive <- round(covFun(lapply(lapply(data,`[[`,"naive"),`[[`,"coverage"),val_naive_sample_coverage[[alg]]),3)[1,]
  cov_wtd <- round(covFun(lapply(lapply(data,`[[`,"wtd"),`[[`,"coverage"),val_label_results$coverage[[alg]]),3)
  # ambiguity
  amb <- ambFun(lapply(data,`[[`,"wtd"),alg=alg)
  # classification performance
  class_measures_naive_pct <- perfBSFun_pct(lapply(data,`[[`,"naive"),data)
  class_measures_wtd <- perfWTDBSFun(lapply(data,`[[`,"wtd"))
  # survival estimation
  data_obs <- lapply(lapply(data,`[[`,"naive"),`[[`,"obs")
  km_naive_pct <- medWTDBSFun(lapply(lapply(data, `[[`,"naive"),`[[`,"km_pred"),data_obs)
  km_wtd <- medWTDBSFun(lapply(lapply(data,`[[`,"wtd"),`[[`,"km_pred"),data_obs)
  out <- list(class_coverage_naive=cov_naive, class_coverage_wtd=cov_wtd,
              ambiguity=amb, class_measures_naive_pct=class_measures_naive_pct,
              class_measures_wtd=class_measures_wtd, km_naive_pct=km_naive_pct,
              km_wtd=km_wtd)
  return(out)
}
bs_measures_list <- map2(bs_out_list, alg_names, ~bsMeasFun(.x,.y))
rm(bs_out_list)

savepath <- c("~/save_path/")

saveRDS(dev_naive_results, paste0(savepath, "dev_naive_fit_algs", ".RDS"))
saveRDS(dev_wtd_results, paste0(savepath, "dev_wtd_fit_algs", ".RDS"))
saveRDS(data_dev_thlds, paste0(savepath, "thlds", ".RDS"))
saveRDS(val_label_results, paste0(savepath, "val_wtd_cov_amb_in_sample", ".RDS"))
saveRDS(val_naive_sample_coverage, paste0(savepath, "naive_sample_coverage", ".RDS"))
saveRDS(val_naive_pred_in_sample, paste0(savepath, "naive_class_perf_in_sample", ".RDS"))
saveRDS(sample_km_obs, paste0(savepath, "km_obs_in_sample", ".RDS"))
saveRDS(sample_km_pred_naive, paste0(savepath, "km_naive_pred_in_sample", ".RDS"))
saveRDS(bs_measures_list, paste0(savepath, "bs_measures_list", ".RDS"))


