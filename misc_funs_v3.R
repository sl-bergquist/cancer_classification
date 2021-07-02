
splitLABELFun <- function(data,...){
  # threshold estimation
  conf_thlds <- function(phat,Y2,class_alphas){
    m = nrow(phat)
    K = ncol(phat)
    score = rep(0,m)
    for(i in 1:m){
      score[i] = phat[i,Y2[i]]
    }
    # dropped total coverage bits 
    
    ## Class-specific coverage
    class_thlds <- rep(NA,K)
    for (k in 1:3){
      class_thlds[k] <- sort(score[Y2==k])[ceiling(class_alphas[k]*(sum(Y2==k)+1)-1)]
    }
    
    return(list(class_thlds=class_thlds))
    
  } 
  
  thlds_est <- conf_thlds(phat=data[,c("p1","p2","p3")], 
                          Y2=data$y_obs,class_alphas=c(alpha,alpha,alpha))
  
  out <- list(thlds_est=thlds_est)
  return(out)
}

LABELFun <- function(pred, thlds,...){
  
  hstarClassFun <- function(pred, thlds){
    tmp <- pred[,c("p1","p2","p3")]
    t(apply(tmp, 1,
            function(x){as.numeric(x >= thlds)}))
  }
  
  Hstar_classcov <- hstarClassFun(pred, thlds)
  colnames(Hstar_classcov) <- c("g1", "g2", "g3")
  
  # ambiguity
  # identify ambiguous Hstars and how many labels they get
  ambFun <- function(x){case_when(
    rowSums(x[,c("g1", "g2", "g3")])==0 ~ "0",
    rowSums(x[,c("g1", "g2", "g3")])==1 ~ "1",
    rowSums(x[,c("g1", "g2", "g3")])==2 ~ "2",
    TRUE ~ "3")}
  label_data <- Hstar_classcov %>% as.data.frame %>% 
    mutate(amb_tmp=ambFun(.),
           amb_flag=factor(amb_tmp,
                           levels=c("0","1","2","3"), ordered=T)) %>% 
    select(-amb_tmp)
  print("label data checkpt")
  # want number of obs per number of labels groups
  amb_flag <- table(label_data$amb_flag) # extract amg_flag, create contingency table (respects factor levels
  amb_flag <- as.data.frame.matrix(t(amb_flag)) # transpose, convert to df
  colnames(amb_flag) <- c("l0", "l1", "l2", "l3") # 0, one, two , or 3 labels
  
  amb_flag_class <- as.data.frame.matrix(table(pred[,"y_obs"],label_data$amb_flag))
  colnames(amb_flag_class) <- c("l0","l1","l2","l3")
  # coverage
  coverage_class <- sapply(1:3,
                           function(k)mean(Hstar_classcov[k==dat_val[,"Y"],k]))
  
  out <- list(Hstar_classcov=Hstar_classcov,amb_flag=amb_flag,
              amb_flag_class=amb_flag_class, coverage_class=coverage_class)
  return(out)
}

## Classification Performance ##
# class1
binCFun1 <- function(df){
  binC1 <- matrix(NA, 2,2)
  binC1[1,1] <- sum(df[,"obs1"]==1&df[,"V1"]==1) # tp
  binC1[1,2] <- sum(df[,"obs1"]==1&df[,"V1"]==0) # fn
  binC1[2,1] <- sum(df[,"obs1"]==0&df[,"V1"]==1) # fp
  binC1[2,2] <- sum(df[,"obs1"]==0&df[,"V1"]==0) #tn
  rownames(binC1) <- c("True_1", "True_0")
  colnames(binC1) <- c("Pred_1", "Pred_0")
  return(binC1)
}
# class 2
binCFun2 <- function(df){
  binC2 <- matrix(NA, 2,2)
  binC2[1,1] <- sum(df[,"obs2"]==1&df[,"V2"]==1) # tp
  binC2[1,2] <- sum(df[,"obs2"]==1&df[,"V2"]==0) # fn
  binC2[2,1] <- sum(df[,"obs2"]==0&df[,"V2"]==1) # fp
  binC2[2,2] <- sum(df[,"obs2"]==0&df[,"V2"]==0) #tn
  rownames(binC2) <- c("True_2", "True_0")
  colnames(binC2) <- c("Pred_2", "Pred_0")
  return(binC2)
}
# class 3
binCFun3 <- function(df){
  binC3 <- matrix(NA, 2,2)
  binC3[1,1] <- sum(df[,"obs3"]==1&df[,"V3"]==1) # tp
  binC3[1,2] <- sum(df[,"obs3"]==1&df[,"V3"]==0) # fn
  binC3[2,1] <- sum(df[,"obs3"]==0&df[,"V3"]==1) # fp
  binC3[2,2] <- sum(df[,"obs3"]==0&df[,"V3"]==0) #tn
  rownames(binC3) <- c("True_3", "True_0")
  colnames(binC3) <- c("Pred_3", "Pred_0")
  return(binC3)
}
# calculate true positives etc. from binary matrices
quadFun <- function(input_mat){
  tp <- input_mat[1,1]
  fp <- input_mat[2,1]
  fn <- input_mat[1,2]
  tn <- input_mat[2,2]
  out <- c(tp=tp, fp=fp, fn=fn, tn=tn)
  return(out)
}
# recall/sensitivity
macRecFun <- function(mat){
  c1 <- mat[1,1]/sum(mat[1,1], mat[3,1]) # tp / (tp+fn)
  c2 <- mat[1,2]/sum(mat[1,2], mat[3,2])
  c3 <- mat[1,3]/sum(mat[1,3], mat[3,3])
  mac_rec <- (c1+c2+c3)/3
  names(mac_rec) <- c("mac_rec")
  mac_rec<- round(mac_rec,4)
  out <- list(mac_rec=mac_rec, c1=c1,c2=c2, c3=c3 )
  return(out)
}
# specificity
macSpecFun <- function(mat){
  c1 <- mat[4,1]/sum(mat[4,1], mat[2,1]) # tn / (tn+fp)
  c2 <- mat[4,2]/sum(mat[4,2], mat[2,2])
  c3 <- mat[4,3]/sum(mat[4,3], mat[2,3])
  mac_spec <- (c1+c2+c3)/3
  names(mac_spec) <- c("mac_spec")
  mac_spec <- round(mac_spec,4)
  out <- list(mac_spec=mac_spec, c1=c1, c2=c2, c3=c3)
  return(out)
}
# precsion/ppv
macPrecFun <- function(mat){
  c1 <- mat[1,1]/sum(mat[1:2,1]) #tp/(tp+fp)
  c2 <- mat[1,2]/sum(mat[1:2,2])
  c3 <- mat[1,3]/sum(mat[1:2,3])
  mac_prec <- (c1+c2+c3)/3
  names(mac_prec) <- c("mac_prec")
  mac_prec <- round(mac_prec,4)
  out <- list(mac_prec=mac_prec, c1=c1, c2=c2, c3=c3)
  return(out)
}
# average accuracy
avgAccFun <- function(mat){
  c1 <- sum(mat[1,1], mat[4,1])/sum(mat[1:4,1])
  c2 <- sum(mat[1,2], mat[4,2])/sum(mat[1:4,2])
  c3 <- sum(mat[1,3], mat[4,3])/sum(mat[1:4,3])
  avgAcc <- (c1+c2+c3)/3
  names(avgAcc) <- c("avgAcc")
  avgAcc <- round(avgAcc,4)
  out <- list(avgAcc=avgAcc, c1=c1, c2=c2, c3=c3)
  return(out)
  
}
# version to implement with weighted bootstrapped data: percentile
perfWTDBSFun <- function(data){
  indFun <- function(x){
    x <- x %>% mutate(V1=ifelse(y_pred_wtd==1,1,0), V2=ifelse(y_pred_wtd==2,1,0), 
                      V3=ifelse(y_pred_wtd==3,1,0),
                      obs1=ifelse(y_obs==1,1,0), 
                      obs2=ifelse(y_obs==2,1,0), obs3=ifelse(y_obs==3,1,0))
    return(x)
  }
  #data <- bs_out
  cmat <- lapply(lapply(data, `[[`, "label_cols"), indFun)
  # matrices of binary results for each class
  bin1 <- lapply(cmat, binCFun1)
  bin2 <- lapply(cmat, binCFun2)
  bin3 <- lapply(cmat, binCFun3)
  
  # calculate true positives etc. from binary matrices
  quad1 <- lapply(bin1, quadFun)
  quad2 <- lapply(bin2, quadFun)
  quad3 <- lapply(bin3, quadFun)
  
  # 3 lists (classes) for each alg
  quad_list <- vector("list", n_bs) # 500 b/c that is number of bs iterations
  for(i in seq_along(1:n_bs)){
    quad_list[[i]] <- cbind(quad1[[i]], quad2[[i]], quad3[[i]]) 
  }
  # recall/sensitivity
  macRecRes <- lapply(quad_list, macRecFun)
  
  # specificity
  macSpecRes <- lapply(quad_list, macSpecFun)
  
  # precsion/ppv
  macPrecRes <- lapply(quad_list, macPrecFun)
  
  # average accuracy
  avgAccRes <- lapply(quad_list, avgAccFun)
  
  ### then need to take upper and lower bounds
  # combine into a data frame across all resamples
  
  tp_list <- lapply(quad_list, `[`,1,)
  fp_list <- lapply(quad_list, `[`,2,)
  fn_list <- lapply(quad_list, `[`,3,)
  tn_list <- lapply(quad_list, `[`,4,)
  
  # tp
  df_tp <- do.call(rbind.data.frame, tp_list)
  colnames(df_tp) <- c("c1","c2","c3")
  tp <- df_tp %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # fp
  df_fp <- do.call(rbind.data.frame, fp_list)
  colnames(df_fp) <- c("c1","c2","c3")
  fp <- df_fp %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # fn
  df_fn <- do.call(rbind.data.frame, fn_list)
  colnames(df_fn) <- c("c1","c2","c3")
  fn <- df_fn %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # tn
  df_tn <- do.call(rbind.data.frame, tn_list)
  colnames(df_tn) <- c("c1","c2","c3")
  tn <- df_tn %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  
  
  df_sens <- do.call(rbind.data.frame, macRecRes)
  
  sens <- df_sens %>% 
    summarise(all=mean(mac_rec, na.rm=T),
          lb_all=quantile(mac_rec, .025, na.rm=T),
           ub_all=quantile(mac_rec, .975, na.rm=T),
          g1=mean(c1, na.rm=T),
           lb_g1=quantile(c1, .025, na.rm=T),
           ub_g1=quantile(c1,.975,na.rm=T),
          g2=mean(c2, na.rm=T),
           lb_g2=quantile(c2,.025,na.rm=T),
           ub_g2=quantile(c2,.975,na.rm=T),
          g3=mean(c3, na.rm=T),
           lb_g3=quantile(c3,.025,na.rm=T),
           ub_g3=quantile(c3,.975,na.rm=T))  
  
  
  
  df_spec <- do.call(rbind.data.frame, macSpecRes)
  
  spec <- df_spec %>% 
    summarise(all=mean(mac_spec, na.rm=T),
              lb_all=quantile(mac_spec, .025, na.rm=T),
              ub_all=quantile(mac_spec, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
   
  
  df_ppv <- do.call(rbind.data.frame, macPrecRes)
  ppv <- df_ppv %>% 
    summarise(all=mean(mac_prec, na.rm=T),
              lb_all=quantile(mac_prec, .025, na.rm=T),
              ub_all=quantile(mac_prec, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  df_acc <- do.call(rbind.data.frame, avgAccRes)
  acc <- df_acc %>% 
    summarise(all=mean(avgAcc, na.rm=T),
              lb_all=quantile(avgAcc, .025, na.rm=T),
              ub_all=quantile(avgAcc, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  class_measures <- list(sens=sens, spec=spec, 
                         ppv=ppv, acc=acc, tp=tp, tn=tn,
                         fp=fp, fn=fn)
  return(class_measures)
}

# version to implement w/ naive bs data: percentile
perfBSFun_pct <- function(data){
  indFun <- function(x){
    x <- x %>% mutate(V1=ifelse(y_pred==1,1,0), V2=ifelse(y_pred==2,1,0), 
                      V3=ifelse(y_pred==3,1,0),
                      obs1=ifelse(y_obs==1,1,0), 
                      obs2=ifelse(y_obs==2,1,0), obs3=ifelse(y_obs==3,1,0))
    return(x)
  }
  #data <- bs_out
  cmat <- lapply(lapply(data, `[[`, "label_cols"), indFun)
  # matrices of binary results for each class
  bin1 <- lapply(cmat, binCFun1)
  bin2 <- lapply(cmat, binCFun2)
  bin3 <- lapply(cmat, binCFun3)
  
  # calculate true positives etc. from binary matrices
  quad1 <- lapply(bin1, quadFun)
  quad2 <- lapply(bin2, quadFun)
  quad3 <- lapply(bin3, quadFun)
  
  # 3 lists (classes) for each alg
  quad_list <- vector("list", n_bs) # 500 b/c that is number of bs iterations
  for(i in seq_along(1:n_bs)){
    quad_list[[i]] <- cbind(quad1[[i]], quad2[[i]], quad3[[i]]) 
  }
  # recall/sensitivity
  macRecRes <- lapply(quad_list, macRecFun)
  
  # specificity
  macSpecRes <- lapply(quad_list, macSpecFun)
  
  # precsion/ppv
  macPrecRes <- lapply(quad_list, macPrecFun)
  
  # average accuracy
  avgAccRes <- lapply(quad_list, avgAccFun)
  
  ### then need to take upper and lower bounds
  # combine into a data frame across all resamples
  
  tp_list <- lapply(quad_list, `[`,1,)
  fp_list <- lapply(quad_list, `[`,2,)
  fn_list <- lapply(quad_list, `[`,3,)
  tn_list <- lapply(quad_list, `[`,4,)
  
  # tp
  df_tp <- do.call(rbind.data.frame, tp_list)
  colnames(df_tp) <- c("c1","c2","c3")
  tp <- df_tp %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # fp
  df_fp <- do.call(rbind.data.frame, fp_list)
  colnames(df_fp) <- c("c1","c2","c3")
  fp <- df_fp %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # fn
  df_fn <- do.call(rbind.data.frame, fn_list)
  colnames(df_fn) <- c("c1","c2","c3")
  fn <- df_fn %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  # tn
  df_tn <- do.call(rbind.data.frame, tn_list)
  colnames(df_tn) <- c("c1","c2","c3")
  tn <- df_tn %>% summarise(g1=mean(c1, na.rm=T),
                            lb_g1=quantile(c1, .025, na.rm=T),
                            ub_g1=quantile(c1, .975, na.rm=T),
                            g2=mean(c2, na.rm=T),
                            lb_g2=quantile(c2, .025, na.rm=T),
                            ub_g2=quantile(c2, .975, na.rm=T),
                            g3=mean(c3, na.rm=T),
                            lb_g3=quantile(c3, .025, na.rm=T),
                            ub_g3=quantile(c3, .975, na.rm=T))
  
  
  df_sens <- do.call(rbind.data.frame, macRecRes)
  
  sens <- df_sens %>% 
    summarise(all=mean(mac_rec, na.rm=T),
              lb_all=quantile(mac_rec, .025, na.rm=T),
              ub_all=quantile(mac_rec, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  
  
  df_spec <- do.call(rbind.data.frame, macSpecRes)
  
  spec <- df_spec %>% 
    summarise(all=mean(mac_spec, na.rm=T),
              lb_all=quantile(mac_spec, .025, na.rm=T),
              ub_all=quantile(mac_spec, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  
  df_ppv <- do.call(rbind.data.frame, macPrecRes)
  ppv <- df_ppv %>% 
    summarise(all=mean(mac_prec, na.rm=T),
              lb_all=quantile(mac_prec, .025, na.rm=T),
              ub_all=quantile(mac_prec, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  df_acc <- do.call(rbind.data.frame, avgAccRes)
  acc <- df_acc %>% 
    summarise(all=mean(avgAcc, na.rm=T),
              lb_all=quantile(avgAcc, .025, na.rm=T),
              ub_all=quantile(avgAcc, .975, na.rm=T),
              g1=mean(c1, na.rm=T),
              lb_g1=quantile(c1, .025, na.rm=T),
              ub_g1=quantile(c1,.975,na.rm=T),
              g2=mean(c2, na.rm=T),
              lb_g2=quantile(c2,.025,na.rm=T),
              ub_g2=quantile(c2,.975,na.rm=T),
              g3=mean(c3, na.rm=T),
              lb_g3=quantile(c3,.025,na.rm=T),
              ub_g3=quantile(c3,.975,na.rm=T))  
  
  class_measures <- list(sens=sens, spec=spec, 
                         ppv=ppv, acc=acc, tp=tp, tn=tn,
                         fp=fp, fn=fn)
  return(class_measures)
}

# version to implement with single sample/iteration
perfFun <- function(data){
  indFun <- function(x){
    x <- x %>% mutate(V1=ifelse(y_pred==1,1,0), V2=ifelse(y_pred==2,1,0), 
                      V3=ifelse(y_pred==3,1,0),
                      obs1=ifelse(y_obs==1,1,0), 
                      obs2=ifelse(y_obs==2,1,0), obs3=ifelse(y_obs==3,1,0))
    return(x)
  }
  
  cmat <- indFun(data)
  # matrices of binary results for each class
 
  bin1 <- binCFun1(cmat)
  bin2 <- binCFun2(cmat)
  bin3 <- binCFun3(cmat)
  
  # calculate true positives etc. from binary matrices

  quad1 <- quadFun(bin1)
  quad2 <- quadFun(bin2)
  quad3 <- quadFun(bin3)
  
  # 3 lists (classes) for each alg
  quad_list <- cbind(quad1, quad2, quad3)
  colnames(quad_list) <- c("g1","g2","g3")
  tp <- quad_list[1,]
  fp <- quad_list[2,]
  fn <- quad_list[3,]
  tn <- quad_list[4,]
  
  # recall/sensitivity
  macRecRes <- macRecFun(quad_list)
  
  # specificity
  macSpecRes <- macSpecFun(quad_list)
  
  # precsion/ppv
  macPrecRes <- macPrecFun(quad_list)
  
  # average accuracy
  avgAccRes <- avgAccFun(quad_list)
  
 
  class_measures <- list(sens=macRecRes, spec=macSpecRes, 
                         ppv=macPrecRes, acc=avgAccRes,
                         tp=tp, tn=tn, fp=fp, fn=fn)
  return(class_measures)
}


# calculate naive standard practice bias measures (bs versions within medWTDBSFun)
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
  return(pct_bias)
}
# version to implement with weighted bootstrap: percentile 
medWTDBSFun <- function(data, data_obs){
  # what strata used in sample?
  names3 <- c("g1","g2","g3")
  
  # check bootstrap iterations -- see if there are fewer than 3
  # this is lazy and doens't distinguish btw classes, but funciton will stop
  # if it's not 1/3 here
  bschkFun <- function(data){
    tmp <- do.call(rbind.data.frame,lapply(lapply(data, `[[`,"n"), length))
    colnames(tmp) <- "count"
    as.list(which(tmp$count<3))
  }
  problem_reps <- bschkFun(data)
  
  # time points
  timeFun <- function(sample,time){summary(sample, times=time, extend=T)$surv}
  # obs version for calculating bias
  obs_30d <- lapply(data_obs,timeFun, time=30)
  obs_90d <- lapply(data_obs,timeFun, time=90)
  obs_365d <- lapply(data_obs, timeFun, time=365)
  obs_50 <- lapply(lapply(data_obs, quantile, probs=.5, conf.int=F),t)
  
  obs_30d <- do.call(rbind.data.frame, obs_30d)
  obs_50 <- do.call(rbind.data.frame, obs_50)
  obs_90d <- do.call(rbind.data.frame, obs_90d)
  obs_365d <- do.call(rbind.data.frame, obs_365d)
  
  tmp_30d <- lapply(data, timeFun, time=30)
  tmp_90d <- lapply(data, timeFun, time=90)
  tmp_365d <- lapply(data, timeFun, time=365)
 
  # median survival time
  #data <- lapply(bs_out,`[[`,"km_pred")
  tmp_50 <- lapply(lapply(data, quantile, probs=.5, conf.int=F),t)
  # find bs reps that have less than 3 classes, fill in empty class with NA
  if (length(problem_reps>0)){
    for(i in problem_reps){
      tmp_50[[i]][3] <- NA
      tmp_30d[[i]][3] <- NA
      tmp_90d[[i]][3] <- NA
      tmp_365d[[i]][3] <- NA
      # tmp_50[[i]] <- c(tmp_50[[i]][1],tmp_50[[i]][3],tmp_50[[i]][2])
    }
  } else(print("no problem reps"))
  tmp_50 <- do.call(rbind.data.frame, tmp_50)
  tmp_30d <- do.call(rbind.data.frame, tmp_30d)
  tmp_90d <- do.call(rbind.data.frame, tmp_90d)
  tmp_365d <- do.call(rbind.data.frame, tmp_365d)
  
  colnames(tmp_30d) <- names3; colnames(tmp_90d) <- names3
  colnames(tmp_365d) <- names3; colnames(tmp_50) <- names3
  colnames(obs_30d) <- names3; colnames(obs_90d) <- names3
  colnames(obs_365d) <- names3; colnames(obs_50) <- names3
  
  dfFun <- function(df){
    tmp <- c(mean_g1=mean(df[["g1"]],na.rm=T),
             lb_g1=quantile(df[["g1"]],.025,na.rm=T),
             ub_g1=quantile(df[["g1"]],.975,na.rm=T),
             mean_g2=mean(df[["g2"]],na.rm=T),
             lb_g2=quantile(df[["g2"]],.025,na.rm=T),
             ub_g2=quantile(df[["g2"]],.975,na.rm=T),
             mean_g3=mean(df[["g3"]],na.rm=T),
             lb_g3=quantile(df[["g3"]],.025,na.rm=T),
             ub_g3=quantile(df[["g3"]],.975,na.rm=T))
    tmp <- as.data.frame(t(tmp))
    colnames(tmp) <- c("g1","lb_g1","ub_g1","g2","lb_g2","ub_g2","g3","lb_g3","ub_g3")
                              
    return(tmp)
  
  }
  df_30d <- dfFun(tmp_30d)
  df_90d <- dfFun(tmp_90d)
  df_365d <- dfFun(tmp_365d)
  df_50 <- dfFun(tmp_50)
  print("dfFun done")
  
  print(dim(tmp_30d));print(head(tmp_30d)); print("dim 30d")
  print(dim(obs_30d));print(head(obs_30d)); print("obs 30d")
  biasFun <- function(pred, obs){
    tmp1 <- pred["g1"] - obs["g1"]
    tmp2 <- pred["g2"] - obs["g2"]
    tmp3 <- pred["g3"] - obs["g3"]
    tmp_df <- cbind(tmp1, tmp2, tmp3)
    tmp_df <- tmp_df %>% summarise(mean_g1=mean(g1, na.rm=T),
                                   lb_g1=quantile(g1,.025,na.rm=T),
                                   ub_g1=quantile(g1,.975,na.rm=T),
                                   mean_g2=mean(g2, na.rm=T),
                                   lb_g2=quantile(g2,.025,na.rm=T),
                                   ub_g2=quantile(g2,.975, na.rm=T),
                                   mean_g3=mean(g3, na.rm=T),
                                   lb_g3=quantile(g3,.025,na.rm=T),
                                   ub_g3=quantile(g3,.975,na.rm=T)) 
    return(tmp_df)
  }
  bias_30d <- biasFun(tmp_30d, obs_30d)
  bias_90d <- biasFun(tmp_90d, obs_90d)
  bias_365d <- biasFun(tmp_365d, obs_365d)
  bias_50 <- biasFun(tmp_50, obs_50)
  
  print("biasFun done")
  
  biasPctFun <- function(pred, obs){
    tmp1 <- 100*((pred["g1"] - obs["g1"])/obs["g1"])
    tmp2 <- 100*((pred["g2"] - obs["g2"])/obs["g2"])
    tmp3 <- 100*((pred["g3"] - obs["g3"])/obs["g3"])
    tmp_df <- cbind(tmp1, tmp2, tmp3)
    tmp_df <- tmp_df %>% summarise(mean_g1=mean(g1, na.rm=T),
                                   lb_g1=quantile(g1,.025,na.rm=T),
                                   ub_g1=quantile(g1,.975,na.rm=T),
                                   mean_g2=mean(g2, na.rm=T),
                                   lb_g2=quantile(g2,.025,na.rm=T),
                                   ub_g2=quantile(g2,.975, na.rm=T),
                                   mean_g3=mean(g3, na.rm=T),
                                   lb_g3=quantile(g3,.025,na.rm=T),
                                   ub_g3=quantile(g3,.975,na.rm=T)) 
    return(tmp_df)
  }
  bias_pct_30d <- biasPctFun(tmp_30d, obs_30d)
  bias_pct_90d <- biasPctFun(tmp_90d, obs_90d)
  bias_pct_365d <- biasPctFun(tmp_365d, obs_365d)
  bias_pct_50 <- biasPctFun(tmp_50, obs_50)
  
  out <- list(d30=df_30d, bias_30=bias_30d, bias_pct_30=bias_pct_30d,
    d90=df_90d,bias_90=bias_90d, bias_pct_90=bias_pct_90d,
              d365=df_365d, bias_365=bias_365d, bias_pct_365=bias_pct_365d,
              median=df_50, bias_med=bias_50, bias_pct_med=bias_pct_50)  
  return(out)
  
}

# version to implement with naive bootstrap: percentile

# coverage
covFun <- function(data, sample_cov){
  # combine into a data frame across all resamples
  df <- do.call(rbind.data.frame, data)
  colnames(df) <- c("g1", "g2", "g3")
  # compute difference between sample cov and bs cov for each class
  df <- df %>% transmute(g1=g1-sample_cov[1],
                         g2=g2-sample_cov[2],
                         g3=g3-sample_cov[3])
  cov <- df %>% 
    summarise(mean_g1=sample_cov[1],
              lb_g1=sample_cov[1]-quantile(g1, .975, na.rm=T),
              ub_g1=sample_cov[1]-quantile(g1, .025, na.rm=T),
              mean_g2=sample_cov[2],
              lb_g2=sample_cov[2]-quantile(g2, .975, na.rm=T),
              ub_g2=sample_cov[2]-quantile(g2, .025, na.rm=T),
              mean_g3=sample_cov[3],
              lb_g3=sample_cov[3]-quantile(g3, .975, na.rm=T),
              ub_g3=sample_cov[3]-quantile(g3, .025, na.rm=T)) 
  return(cov[1,])
}

# ambiguity
ambFun <- function(data){
  #data <- lapply(bs_out, `[[`, "wtd")
  sample_amb <- as.matrix(label_results$amb_flag)
  label_cols <- lapply(data, `[[`, "label_cols")
  #  extract amg_flag, create contingency table (respects factor levels)
  amb_flag <- lapply(lapply(label_cols, `[[`, "amb_flag"), table) 
  amb_flag <- lapply(lapply(amb_flag, t), as.data.frame.matrix) # transpose, convert to df
  df <- do.call(rbind.data.frame, amb_flag)
  colnames(df) <- c("l0", "l1", "l2", "l3") # 0, one, two , or 3 labels
  df <- df %>% transmute(l0=l0-sample_amb[1],
                      l1=l1-sample_amb[2],
                      l2=l2-sample_amb[3],
                      l3=l3-sample_amb[4])
  
  amb <- df %>% summarise(mean_l0=sample_amb[1],
                          lb_l0=sample_amb[1]-quantile(l0, .975),
                          ub_l0=sample_amb[1]-quantile(l0, .025),
                          mean_l1=sample_amb[2],
                          lb_l1=sample_amb[2]-quantile(l1, .975),
                          ub_l1=sample_amb[2]-quantile(l1, .025),
                          mean_l2=sample_amb[3],
                          lb_l2=sample_amb[3]-quantile(l2, .975),
                          ub_l2=sample_amb[3]-quantile(l2, .025),
                          mean_l3=sample_amb[4],
                          lb_l3=sample_amb[4]-quantile(l3, .975),
                          ub_l3=sample_amb[4]-quantile(l3, .025))
  
  # so this give us total 
  # we also want to capture ambiguity by true class
  amb_flag <- lapply(label_cols, `[`, c("y_obs","amb_flag"))
  # so I just need to do this in the sample version. annoying. 
  amb_table <- lapply(amb_flag, table)
  sample_amb_class <- label_results$amb_flag_class
  # now how to combine across reps? make them single class dfs?
  classFun <- function(class){
    tmp <- lapply(amb_table, function(x) x[class,])
    df <- do.call(rbind.data.frame, tmp)
    colnames(df) <- c("l0", "l1", "l2", "l3")
    df <- df %>% transmute(l0=l0-sample_amb_class[class,1],
                           l1=l1-sample_amb_class[class,2],
                           l2=l2-sample_amb_class[class,3],
                           l3=l3-sample_amb_class[class,4])
    amb_class <- df %>% summarise(mean_l0=sample_amb_class[class,1],
                                  lb_l0=sample_amb_class[class,1]-quantile(l0,.975),
                                  ub_l0=sample_amb_class[class,1]-quantile(l0,.025),
                                  mean_l1=sample_amb_class[class,2],
                                  lb_l1=sample_amb_class[class,2]-quantile(l1,.975),
                                  ub_l1=sample_amb_class[class,2]-quantile(l1,.025),
                                  mean_l2=sample_amb_class[class,3],
                                  lb_l2=sample_amb_class[class,3]-quantile(l2,.975),
                                  ub_l2=sample_amb_class[class,3]-quantile(l2,.025),
                                  mean_l3=sample_amb_class[class,4],
                                  lb_l3=sample_amb_class[class,4]-quantile(l3,.975),
                                  ub_l3=sample_amb_class[class,4]-quantile(l3,.025))
    return(amb_class)
  }
  amb_1 <- classFun(1); amb_2 <- classFun(2); amb_3 <- classFun(3)
  out <- list(amb_all=amb, amb_g1=amb_1, amb_g2=amb_2,amb_g3=amb_3)
  return(out)
}

