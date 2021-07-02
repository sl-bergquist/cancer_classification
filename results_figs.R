##################################################
# CancerCLAS Multiclass Prediction
# Simulation and Data Results
# Figures
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


# colourblind palette w/ black
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#AED2E6")
# shapes
shapePalette <- c(15,19,17,23,8,25,4)


sim_path <- "insert_save_path_string"
setting_label <- c("Scenario 1: Accurate, Certain", "Scenario 2: Accurate, Uncertain",
                   "Scenario 3: Inaccurate, Uncertain")

#### LABEL SET AMBIGUITY ####
amb_all <- read.table(paste0(sim_path, "amb_all", ".txt"), sep=",", header=T)
amb_all_pct <- amb_all %>% mutate_at(vars(-setting), funs(./1000)) %>% mutate_at(vars(-setting), ~round(.,2))

amb_g1 <- read.table(paste0(sim_path, "amb_g1", ".txt"), sep=",", header=T)
amb_g2 <- read.table(paste0(sim_path, "amb_g2", ".txt"), sep=",", header=T)
amb_g3 <- read.table(paste0(sim_path, "amb_g3", ".txt"), sep=",", header=T)
# by class
n_g1 <- (amb_g1$l0 + amb_g1$l1 + amb_g1$l2 + amb_g1$l3)[1]
n_g2 <- (amb_g2$l0 + amb_g2$l1 + amb_g2$l2 + amb_g2$l3)[1]
n_g3 <- (amb_g3$l0 + amb_g3$l1 + amb_g3$l2 + amb_g3$l3)[1]
amb_g1_pct <- amb_g1 %>% mutate_at(vars(-setting), funs(./n_g1)) %>% mutate_at(vars(-setting), ~round(.,2))
amb_g2_pct <- amb_g2 %>% mutate_at(vars(-setting), funs(./n_g2)) %>% mutate_at(vars(-setting), ~round(.,2))
amb_g3_pct <- amb_g3 %>% mutate_at(vars(-setting), funs(./n_g3)) %>% mutate_at(vars(-setting), ~round(.,2))

# make grp var and put these all into a df
amb_all_pct <- amb_all_pct %>% mutate(grp="All", setting=setting_label)
amb_g1_pct <- amb_g1_pct %>% mutate(grp="Class 1", setting=setting_label)
amb_g2_pct <- amb_g2_pct %>% mutate(grp="Class 2", setting=setting_label)
amb_g3_pct <- amb_g3_pct %>% mutate(grp="Class 3", setting=setting_label)

amb_df <- rbind.data.frame(amb_all_pct, amb_g1_pct, amb_g2_pct, amb_g3_pct)
tmp <- amb_df %>% gather("label","percent",-setting, -grp)
tmp_lb <- tmp %>% filter(label %in% c("lb_l0", "lb_l1","lb_l2","lb_l3")) %>% 
  rename(lb=percent) %>% mutate(label=str_sub(label,4,5))
tmp_ub <- tmp %>% filter(label %in% c("ub_l0","ub_l1","ub_l2","ub_l3")) %>%
  rename(ub=percent) %>% mutate(label=str_sub(label,4,5))

tmp <- tmp %>% filter(label %in% c("l0","l1","l2","l3"))
tmp <- tmp %>% left_join(tmp_lb) %>% left_join(tmp_ub)
tmp$label <- recode(tmp$label, "l0"="0","l1"="1","l2"="2","l3"="3")


ggplot(tmp, aes(fill=label, y=percent*100,ymin=lb*100, ymax=ub*100, x=grp)) + 
  geom_bar(position="dodge", stat="identity",width=.9) +
  geom_errorbar(width=.5,position=position_dodge(.9)) +
  scale_fill_manual(values=okabe) +
  scale_y_continuous(expand=c(0,.01))+
 coord_cartesian(ylim=c(0,105))+
  facet_wrap(~setting) +
  theme(panel.border=element_rect(colour="black",fill=NA),
        legend.position="bottom",
        axis.ticks=element_blank(),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(vjust=0),
        strip.text=element_text(size=15),
        panel.background = element_rect(fill="grey92"),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        plot.title=element_text(size=15),
        axis.title.y=element_text(size=15),
        axis.title.x=element_blank())+
  labs(fill="Number of Labels in Set") +
  xlab("") +
  ylab("Percent of Observations") + NULL
  #ggtitle("Simulation Study: Label Set Ambiguity by True Class")

# change to 800 dpi for final versions
# Figure 2: Simulation study label ambiguity: Share of sample by number of assigned 
  # labels in label set.
ggsave(paste0(sim_path,"figures/", "sim_label_amb.tiff"), width=12, height=5.5, units=c("in"), dpi=500)
rm(list=ls(pattern="amb"))
rm(list=ls(pattern="tmp"))

#### data label ambiguity ####
data_path <- "W:/DATA/cms_ocm/Machine_Learning/lung/analysis/bergquist/2012_2013/output/emp/032621/"
amb_all <- read.table(paste0(data_path, "bs_wtd_amb_total_pct", ".txt"), sep=",", header=T)
amb_g1 <- read.table(paste0(data_path, "bs_wtd_amb_g1_pct", ".txt"), sep=",", header=T)
amb_g2 <- read.table(paste0(data_path, "bs_wtd_amb_g2_pct", ".txt"), sep=",", header=T)
amb_g3 <- read.table(paste0(data_path, "bs_wtd_amb_g3_pct", ".txt"), sep=",", header=T)

# make grp var and put these all into a df
amb_all <- amb_all %>% mutate(grp="All")
amb_g1 <- amb_g1 %>% mutate(grp="Stage I/II")
amb_g2 <- amb_g2 %>% mutate(grp="Stage III")
amb_g3 <- amb_g3 %>% mutate(grp="Stage IV")

amb_df <- rbind.data.frame(amb_all, amb_g1, amb_g2, amb_g3)
tmp <- amb_df %>% gather("label","percent",-alg, -grp)
tmp_lb <- tmp %>% filter(label %in% c("lb_l0", "lb_l1","lb_l2","lb_l3")) %>% 
  rename(lb=percent) %>% mutate(label=str_sub(label,4,5))
tmp_ub <- tmp %>% filter(label %in% c("ub_l0","ub_l1","ub_l2","ub_l3")) %>%
  rename(ub=percent) %>% mutate(label=str_sub(label,4,5))

tmp <- tmp %>% mutate(label=str_sub(label,6,7)) %>% filter(label %in% c("l0","l1","l2","l3"))
tmp <- tmp %>% left_join(tmp_lb) %>% left_join(tmp_ub)
tmp$label <- recode(tmp$label, "l0"="0","l1"="1","l2"="2","l3"="3")
tmp$alg <- recode(tmp$alg, "enet"="Elastic Net Regression", "gam"="Generalized Additive Regression",
                  "lasso"="Lasso Regression","mn"="Multinomial Logistic Regression",
                  "rf"="Random Forests", "ridge"="Ridge Regression", "xgb"="Gradient Boosting")

tmp <- tmp %>% filter(label!=0)
## selected algs
tmp <- tmp %>% filter(alg %in% c("Multinomial Logistic Regression", "Random Forests"))

ggplot(tmp, aes(fill=label, y=percent*100,ymin=lb*100,ymax=ub*100, x=grp)) + 
  geom_bar(position="dodge", stat="identity",width=.9) +
  coord_cartesian(ylim=c(0,100))+
  geom_errorbar(width=.5, position=position_dodge(width=,.9)) +
  scale_fill_manual(values=okabe[2:4]) +
  scale_y_continuous(expand=c(0,.01))+
  facet_wrap(~alg) +
  theme(legend.position="bottom",
        legend.direction = "horizontal",
        panel.background = element_rect(fill="grey92"),
        panel.border=element_rect(colour="black",fill=NA),
        panel.grid.major.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        axis.ticks=element_blank(),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.title=element_text(size=15),
        legend.text=element_text(size=15),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(vjust=0),
        strip.text=element_text(size=15),
        plot.title=element_text(size=15),
        axis.title.y=element_text(size=15)) +
  labs(fill="Number of Labels in Set") +
  xlab("") +
  ylab("Percent of Observations") + NULL  
#ggtitle("Data Analysis: Label Set Ambiguity by True Class (Selected Algorithms)") +

# Figure 4: Data analysis label ambiguity: Share of sample by number of assigned 
  # labels in label set
ggsave(paste0(data_path,"figures/", "data_label_amb_selected.tiff"), width=12, height=5.5, units=c("in"), dpi=500)


setting_label2 <- c("Scenario 1:\nAccurate, Certain", "Scenario 2:\nAccurate, Uncertain",
                   "Scenario 3:\nInaccurate, Uncertain")

#### SIM MEDIAN SURVIVAL BIAS DAYS ####

# facet grid
surv <- read.table(paste0(sim_path, "median_bias_mean", ".txt"), sep=",", header=T)
surv_lb <- read.table(paste0(sim_path, "median_bias_lb", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class", "lb",-setting)
surv_ub <- read.table(paste0(sim_path, "median_bias_ub", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class","ub",-setting)
surv_df <- surv %>% gather("class","value",-setting)  %>% 
  left_join(surv_lb) %>% left_join(surv_ub)

surv_df <- surv_df %>% mutate(method=case_when(
  setting %in% c("naive_bs_set_1","naive_bs_set_2","naive_bs_set_3") ~ "Naive Bootstrap",
  setting %in% c("wtd_bs_set_1","wtd_bs_set_2","wtd_bs_set_3") ~ "Weighted Bootstrap",
  TRUE ~ "Naive Standard Practice")) %>% 
  mutate(setting=case_when(
    setting %in% c("naive_sample_set_1","naive_bs_set_1","wtd_bs_set_1") ~ setting_label2[1],
    setting %in% c("naive_sample_set_2", "naive_bs_set_2","wtd_bs_set_2") ~ setting_label2[2],
    TRUE ~ setting_label2[3])) %>% 
  mutate(method=factor(method, ordered=T, 
                       levels=c("Naive Standard Practice",
                                "Naive Bootstrap",
                                "Weighted Bootstrap"))) %>%
  mutate(class=case_when(class=="g1"~"Class 1",class=="g2"~"Class 2",TRUE~"Class 3"))
  
# suppress really tight CIs: days -- suppress below 31 days
tmp <- surv_df %>% mutate(lb=ifelse(ub-lb<31,NA,lb),
                              ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=method, y=value,shape=method,colour=setting)) +
  geom_errorbar(aes(x=method, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +
 # coord_cartesian(ylim=c(-100,150))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  facet_grid(class ~ setting) +
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("Median Survival: Bias (Days)")
# Figure 3: Simulation study median survival days bias
ggsave(paste0(sim_path,"figures/", "sim_med_bias.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

### BIAS FOR 90 AND 365 PROBABILITY - SIM STUDY ###

#90 day survival probability bias
# facet grid
surv <- read.table(paste0(sim_path, "d90_bias_mean", ".txt"), sep=",", header=T)
surv_lb <- read.table(paste0(sim_path, "d90_bias_lb", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class", "lb",-setting)
surv_ub <- read.table(paste0(sim_path, "d90_bias_ub", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class","ub",-setting)
surv_df <- surv %>% gather("class","value",-setting)  %>% 
  left_join(surv_lb) %>% left_join(surv_ub)

surv_df <- surv_df %>% mutate(method=case_when(
  setting %in% c("naive_bs_set_1","naive_bs_set_2","naive_bs_set_3") ~ "Naive Bootstrap",
  setting %in% c("wtd_bs_set_1","wtd_bs_set_2","wtd_bs_set_3") ~ "Weighted Bootstrap",
  TRUE ~ "Naive Standard Practice")) %>% 
  mutate(setting=case_when(
    setting %in% c("naive_set_1","naive_bs_set_1","wtd_bs_set_1") ~ setting_label2[1],
    setting %in% c("naive_set_2", "naive_bs_set_2","wtd_bs_set_2") ~ setting_label2[2],
    TRUE ~ setting_label2[3])) %>% 
  mutate(method=factor(method, ordered=T, 
                       levels=c("Naive Standard Practice",
                                "Naive Bootstrap",
                                "Weighted Bootstrap"))) %>%
  mutate(class=case_when(class=="g1"~"Class 1",class=="g2"~"Class 2",TRUE~"Class 3"))

# suppress really tight CIs: probability (bias) so below 0.05
tmp <- surv_df %>% mutate(lb=ifelse(ub-lb<.05,NA,lb),
                          ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=method, y=value,shape=method,colour=setting)) +
  geom_errorbar(aes(x=method, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +
  coord_cartesian(ylim=c(-.26,1.05))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_y_continuous(breaks=c(-0.25, 0.00, 0.25,0.50,0.75,1.00))+
  facet_grid(class ~ setting) +
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("90-day Survival Probability: Bias")
# Figure B2: Simulation study 90-day survival probability bias
ggsave(paste0(sim_path,"figures/", "sim_d90_bias.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

#365 day survival probability bias

# facet grid
surv <- read.table(paste0(sim_path, "d365_bias_mean", ".txt"), sep=",", header=T)
surv_lb <- read.table(paste0(sim_path, "d365_bias_lb", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class", "lb",-setting)
surv_ub <- read.table(paste0(sim_path, "d365_bias_ub", ".txt"), sep=",", header=T) %>%
  select(-contains("sd")) %>% gather("class","ub",-setting)
surv_df <- surv %>% gather("class","value",-setting)  %>% 
  left_join(surv_lb) %>% left_join(surv_ub)

surv_df <- surv_df %>% mutate(method=case_when(
  setting %in% c("naive_bs_set_1","naive_bs_set_2","naive_bs_set_3") ~ "Naive Bootstrap",
  setting %in% c("wtd_bs_set_1","wtd_bs_set_2","wtd_bs_set_3") ~ "Weighted Bootstrap",
  TRUE ~ "Naive Standard Practice")) %>% 
  mutate(setting=case_when(
    setting %in% c("naive_set_1","naive_bs_set_1","wtd_bs_set_1") ~ setting_label2[1],
    setting %in% c("naive_set_2", "naive_bs_set_2","wtd_bs_set_2") ~ setting_label2[2],
    TRUE ~ setting_label2[3])) %>% 
  mutate(method=factor(method, ordered=T, 
                       levels=c("Naive Standard Practice",
                                "Naive Bootstrap",
                                "Weighted Bootstrap"))) %>%
  mutate(class=case_when(class=="g1"~"Class 1",class=="g2"~"Class 2",TRUE~"Class 3"))

# suppress really tight CIs: probabilities, so below 0.05
tmp <- surv_df %>% mutate(lb=ifelse(ub-lb<.05,NA,lb),
                          ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=method, y=value,shape=method,colour=setting)) +
  geom_errorbar(aes(x=method, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +
  # coord_cartesian(ylim=c(-100,150))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  facet_grid(class ~ setting) +
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("365-day Survival Probability: Bias")
# Figure B3: Simulation study 365-day survival probability bias
ggsave(paste0(sim_path,"figures/", "sim_d365_bias.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)



### DATA: MEDIAN SURVIVAL BIAS ###
# median bias --days
surv_bs_wtd <- read.table(paste0(data_path, "bs_wtd_med_bias", ".txt"), sep=",", header=T) 
surv_bs_naive <- read.table(paste0(data_path,"bs_naive_pct_med_bias",".txt"),sep=",",header=T)
surv_naive_sample <- read.table(paste0(data_path,"naive_sample_med_bias",".txt"),sep=",",header=T)

# separate out bounds then gather
surv_bs_wtd_mean <- surv_bs_wtd %>% select(alg, g1=mean_g1, g2=mean_g2, g3=mean_g3) %>%
  gather("class","bias",-alg)
surv_bs_wtd_lb <- surv_bs_wtd %>% select(alg, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class", "lb",-alg)
surv_bs_wtd_ub <- surv_bs_wtd %>% select(alg, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class", "ub",-alg)
surv_bs_wtd <- surv_bs_wtd_mean %>% left_join(surv_bs_wtd_lb) %>% left_join(surv_bs_wtd_ub) %>%
  mutate(method="Weighted Bootstrap")

surv_bs_naive_mean <- surv_bs_naive %>% select(alg, g1=mean_g1, g2=mean_g2, g3=mean_g3) %>% gather("class","bias",-alg)
surv_bs_naive_lb <- surv_bs_naive %>% select(alg, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class", "lb",-alg)
surv_bs_naive_ub <- surv_bs_naive %>% select(alg, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class", "ub",-alg)
surv_bs_naive <- surv_bs_naive_mean %>% left_join(surv_bs_naive_lb) %>% left_join(surv_bs_naive_ub) %>%
  mutate(method="Naive Bootstrap") 

surv_naive_sample <- surv_naive_sample %>% gather("class","bias",-alg) %>%
  mutate(lb=NA, ub=NA) %>%
  mutate(method="Naive Standard Practice")

# rbind
surv_bias_plot <- rbind(surv_bs_wtd, surv_bs_naive, surv_naive_sample)

surv_bias_plot$alg <- recode(surv_bias_plot$alg, "enet"="Elastic Net", "gam"="Generalized Additive",
                  "lasso"="Lasso","mn"="Multinomial Logistic",
                  "rf"="Random Forests", "ridge"="Ridge", "xgb"="Gradient Boosting")
surv_bias_plot$class <- recode(surv_bias_plot$class, "g1"="Stage I/II","g2"="Stage III","g3"="Stage IV")
surv_bias_plot <- surv_bias_plot %>% mutate(method=factor(method, ordered=T, 
                     levels=c("Naive Standard Practice",
                              "Naive Bootstrap",
                              "Weighted Bootstrap")))
# CIs -- dauys, so less than 31
tmp <- surv_bias_plot %>% mutate(lb=ifelse(ub-lb<31,NA,lb),
                                 ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=alg, y=bias,shape=alg,colour=method)) +
  geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  #coord_cartesian(ylim=c(-100,150))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
 #scale_y_continuous(labels = scales::percent)+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("Median Survival: Bias (Days)")
# Figure C6: Data analysis: Median days survival bias
ggsave(paste0(data_path,"figures/", "data_med_bias.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)


### DATA: SURVIVAL PROBABILITY BIAS - 90 and 365 DAYS ###
# survival probability at 90 days (bias) -- pct
biasFun <- function(time){
  surv_bs_wtd <- read.table(paste0(data_path, "bs_wtd_", time, "_bias_pct", ".txt"), sep=",", header=T) 
  surv_bs_naive <- read.table(paste0(data_path,"bs_naive_pct_",time,"_bias_pct",".txt"),sep=",",header=T)
  surv_naive_sample <- read.table(paste0(data_path,"naive_sample_",time, "_bias_pct",".txt"),sep=",",header=T)
  
  
  # separate out bounds then gather
  surv_bs_wtd_mean <- surv_bs_wtd %>% select(alg, g1=mean_g1, g2=mean_g2, g3=mean_g3) %>%
    gather("class","bias",-alg)
  surv_bs_wtd_lb <- surv_bs_wtd %>% select(alg, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class", "lb",-alg)
  surv_bs_wtd_ub <- surv_bs_wtd %>% select(alg, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class", "ub",-alg)
  surv_bs_wtd <- surv_bs_wtd_mean %>% left_join(surv_bs_wtd_lb) %>% left_join(surv_bs_wtd_ub) %>%
    mutate(method="Weighted Bootstrap")
  
  surv_bs_naive_mean <- surv_bs_naive %>% select(alg, g1=mean_g1, g2=mean_g2, g3=mean_g3) %>% gather("class","bias",-alg)
  surv_bs_naive_lb <- surv_bs_naive %>% select(alg, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class", "lb",-alg)
  surv_bs_naive_ub <- surv_bs_naive %>% select(alg, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class", "ub",-alg)
  surv_bs_naive <- surv_bs_naive_mean %>% left_join(surv_bs_naive_lb) %>% left_join(surv_bs_naive_ub) %>%
    mutate(method="Naive Bootstrap")
  
  surv_naive_sample <- surv_naive_sample %>% gather("class","bias",-alg) %>%
    mutate(lb=NA, ub=NA) %>%
    mutate(method="Naive Standard Practice")
  
  # rbind
  surv_bias_plot <- rbind(surv_bs_wtd, surv_bs_naive, surv_naive_sample)
  
  surv_bias_plot$alg <- recode(surv_bias_plot$alg, "enet"="Elastic Net", "gam"="Generalized Additive",
                               "lasso"="Lasso","mn"="Multinomial Logistic",
                               "rf"="Random Forests", "ridge"="Ridge", "xgb"="Gradient Boosting")
  surv_bias_plot$class <- recode(surv_bias_plot$class, "g1"="Stage I/II","g2"="Stage III","g3"="Stage IV")
  surv_bias_plot <- surv_bias_plot %>% mutate(method=factor(method, ordered=T, 
                                                            levels=c("Naive Standard Practice",
                                                                     "Naive Bootstrap",
                                                                     "Weighted Bootstrap")))
   return(surv_bias_plot)
}
surv_plot_d90 <- biasFun("d90")
# suppress tight CIs: <1% (lets do 1% for percent bias)
tmp <- surv_plot_d90 %>% 
  mutate(lb=ifelse(ub-lb<1,NA,lb),
         ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=alg, y=bias,shape=alg,colour=method)) +
  geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(-6,10))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  #scale_y_continuous(labels = scales::percent)+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("90-Day Survival Probability: Percent Bias")
# Figure 6: Data analysis: 90-day survival probability percent bias
ggsave(paste0(data_path,"figures/", "data_d90_bias_pct.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

# 365 days
surv_plot_d365 <- biasFun("d365")
# suppress tight CIs: bootstrap methods
tmp <- surv_plot_d365 %>% 
  mutate(lb=ifelse(ub-lb<1,NA,lb),
         ub=ifelse(is.na(lb),NA,ub))
ggplot(tmp, aes(x=alg, y=bias,shape=alg,colour=method)) +
  geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
 # coord_cartesian(ylim=c(-.2,.2))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("365-Day Survival Probability: Percent Bias")
# Figure C5: Data analysis: 365-day survival probability percent bias
ggsave(paste0(data_path,"figures/", "data_d365_bias_pct.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)


# combo many discrimation metrics figures
### SIMULATION DISCRIMINATION ###
acc <- read.table(paste0(sim_path, "acc", ".txt"), sep=",", header=T) %>%
  mutate(measure="Accuracy")
sens <- read.table(paste0(sim_path, "sens", ".txt"),sep=",",header=T) %>%
  mutate(measure="Sensitivity")
spec <- read.table(paste0(sim_path, "spec", ".txt"),sep=",",header=T) %>%
  mutate(measure="Specificity")
ppv <- read.table(paste0(sim_path, "ppv", ".txt"),sep=",",header=T) %>%
  mutate(measure="PPV")

acc_lb <- acc %>% select(setting, measure, contains("lb")) %>%
  rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
  gather("grp","lb",-setting,-measure)
acc_ub <- acc %>% select(setting, measure, contains("ub")) %>%
  rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
  gather("grp","ub",-setting,-measure)
sens_lb <- sens %>% select(setting, measure, contains("lb")) %>%
  rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
  gather("grp","lb",-setting,-measure)
sens_ub <- sens %>% select(setting, measure, contains("ub")) %>%
  rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
  gather("grp","ub",-setting,-measure)
spec_lb <- spec %>% select(setting, measure, contains("lb")) %>%
  rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
  gather("grp","lb",-setting,-measure)
spec_ub <- spec %>% select(setting, measure, contains("ub")) %>%
  rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
  gather("grp","ub",-setting,-measure)
ppv_lb <- ppv %>% select(setting, measure, contains("lb")) %>%
  rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
  gather("grp","lb",-setting,-measure)
ppv_ub <- ppv %>% select(setting, measure, contains("ub")) %>%
  rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
  gather("grp","ub",-setting,-measure)

acc <- acc %>% select(setting, measure, all, g1, g2, g3) %>% 
  gather("grp","mean",-setting,-measure) %>% left_join(acc_lb) %>% 
  left_join(acc_ub)
sens <- sens %>% select(setting, measure, all, g1, g2, g3) %>% 
  gather("grp","mean",-setting,-measure) %>% left_join(sens_lb) %>% 
  left_join(sens_ub)
spec <- spec %>% select(setting, measure, all, g1, g2, g3) %>% 
  gather("grp","mean",-setting,-measure) %>% left_join(spec_lb) %>% 
  left_join(spec_ub)
ppv <- ppv %>% select(setting, measure, all, g1, g2, g3) %>% 
  gather("grp","mean",-setting,-measure) %>% left_join(ppv_lb) %>% 
  left_join(ppv_ub)

cm_df <- rbind(acc,sens,spec,ppv)

cm_df$method <- recode(cm_df$setting,"naive_xbs_set_1"="Naive Bootstrap",
                       "naive_xbs_set_2"="Naive Bootstrap",
                       "naive_xbs_set_3"="Naive Bootstrap",
                       "naive_set_1"="Naive Standard Practice", "naive_set_2"="Naive Standard Practice",
                       "naive_set_3"="Naive Standard Practice",
                       "wtd_bs_set_1"="Weighted Bootstrap", "wtd_bs_set_2"="Weighted Bootstrap",
                       "wtd_bs_set_3"="Weighted Bootstrap")
cm_df <- cm_df %>% 
  mutate(setting=case_when(
    setting %in% c("naive_set_1","naive_xbs_set_1","wtd_bs_set_1") ~ setting_label2[1],
    setting %in% c("naive_set_2", "naive_xbs_set_2","wtd_bs_set_2") ~ setting_label2[2],
    TRUE ~ setting_label2[3])) %>% 
  mutate(method=factor(method, ordered=T, 
                       levels=c("Naive Standard Practice",
                                "Naive Bootstrap",
                                "Weighted Bootstrap"))) 
cm_df$grp <- recode(cm_df$grp, "g1"="Class 1","g2"="Class 2","g3"="Class 3","all"="Average")


cm_df_avg <- cm_df %>% filter(grp=="Average") # suppress CIs smaller than 5% (probabilities)
cm_df_avg <- cm_df_avg %>% mutate(lb=ifelse(ub-lb<.05,NA,lb),
                                  ub=ifelse(is.na(lb),NA,ub))
ggplot(cm_df_avg, aes(x=method, y=mean*100,shape=method,colour=setting)) +
  geom_errorbar(aes(x=method, ymin=lb*100, ymax=ub*100), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(0,100))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
 # scale_shape_manual(values=shapePalette)+
 # scale_y_continuous(labels = scales::percent)+
  facet_grid(measure ~ setting) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        panel.background = element_rect(fill="grey92"),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("Percent")
# Figure B1: Simulation study classification discrimination: Average accuracy,
  # sensitivity, specificity, PPV
ggsave(paste0(sim_path,"figures/", "sim_cm_avg_pct.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)



### CALIBRATION -- DATA ###
cal_naive_sample <- read.table(paste0(data_path,"val_naive_in_sample" ,"_cal", ".txt"), sep=",", header=T) 
cal_naive_pct <- read.table(paste0(data_path,"bs_naive_pct" ,"_cal", ".txt"), sep=",", header=T) %>%
  select(-lb, -ub) %>% rename(prop_obs=mean_prop_obs)
cal_bs_wtd <- read.table(paste0(data_path,"bs_wtd" ,"_cal", ".txt"), sep=",", header=T) %>%
  select(-lb, -ub) %>% rename(prop_obs=mean_prop_obs)
cal_plot <- rbind.data.frame(cal_naive_sample,cal_naive_pct,cal_bs_wtd)
# only show naive sample predictions -- the other needs fixing, and they are about the same
cal_plot <- cal_plot %>% filter(method!="Naive Bootstrap")
cal_plot$method <- recode(cal_plot$method, "Naive Standard Practice"="Naive")
cal_plot <- cal_plot %>% mutate(class=case_when(
  class=="g1" ~ "Stage I/II",
  class=="g2" ~ "Stage III",
  TRUE ~ "Stage IV"
))
cal_plot$alg <- recode(cal_plot$alg, "enet"="Elastic\nNet", "gam"="Generalized\nAdditive",
                         "lasso"="Lasso","mn"="Multinomial\nLogistic",
                         "rf"="Random\nForests", "ridge"="Ridge", "xgb"="Gradient\nBoosting")


# then group by method: all in same thing
cal_plot %>% ggplot(aes(x=bin, y=prop_obs, colour=method, shape=method)) + 
  geom_point(size=2, alpha=.75) +
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_x_discrete(limits=c(1, 2, 3, 4, 5,6,7,8,9,10),
                   labels=c("10","20","30","40","50",
                            "60","70","80","90 ", "100"), expand=c(0,0))+
  scale_y_continuous(breaks=seq(0,1,.25),
                     labels=c("0", "25", "50", "75", "100")) +
  expand_limits(x=c(0,10))+
  facet_grid(alg ~ class) +
  theme(panel.border=element_rect(colour="black",fill=NA),
        panel.background=element_rect(fill="grey92"),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.position="bottom",
        legend.key =element_blank(),
        legend.text=element_text(size=13),
        legend.title=element_blank(),
        strip.background=element_rect(colour="black",size=.9, fill=okabe[8]),
        strip.text=element_text(size=15),
        axis.text.x=element_text(hjust=1),
        axis.title.x=element_text(vjust=-1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.y = unit(1,"lines"),
        axis.title=element_text(size=15)) +
  geom_vline(xintercept=0, color="grey")+
  labs(x="Predicted Probability Percentile Bin",
       y="Percent Observed Stage") 
# Figure C4: Data analysis: Observed stage by predicted probability
ggsave(paste0(data_path,"figures/", "data_cal.tiff"), width=8.25, height=11, units=c("in"), dpi=500)



### DATA DISCRIMINATION GRID ### 
# maybe gonna need to make this some funcitons -- have 3 measures and 3 methods
cmReadFun <- function(name){
  
  acc <- read.table(paste0(data_path,name ,"_acc", ".txt"), sep=",", header=T) %>%
    mutate(measure="Accuracy")
  sens <- read.table(paste0(data_path,name, "_sens", ".txt"),sep=",",header=T) %>%
    mutate(measure="Sensitivity")
  spec <- read.table(paste0(data_path,name, "_spec", ".txt"),sep=",",header=T) %>%
    mutate(measure="Specificity")
  ppv <- read.table(paste0(data_path,name, "_ppv", ".txt"),sep=",",header=T) %>%
    mutate(measure="PPV")
  
  if(name=="val_naive_in_sample"){
    
    acc <- acc %>% rename(all=avgAcc, g1=c1, g2=c2,g3=c3)
    acc_lb <- acc %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp", "lb",-alg,-measure)
    acc_ub <- acc %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","ub",-alg,-measure)
    
    sens <- sens %>% rename(all=mac_rec,g1=c1,g2=c2,g3=c3)
    sens_lb <- sens %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","lb",-alg,-measure)
    sens_ub <- sens %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","ub",-alg,-measure)
    
    spec <- spec %>% rename(all=mac_spec,g1=c1,g2=c2,g3=c3)
    spec_lb <- spec %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","lb",-alg,-measure)
    spec_ub <- spec %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","ub",-alg,-measure)
    
    ppv <- ppv %>% rename(all=mac_prec,g1=c1,g2=c2,g3=c3)
    ppv_lb <- ppv %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","lb",-alg,-measure)
    ppv_ub <- ppv %>% select(alg, measure) %>% mutate(all=NA, g1=NA, g2=NA, g3=NA) %>%
      gather("grp","ub",-alg,-measure)
    
  }else{
  acc_lb <- acc %>% select(alg, measure, contains("lb")) %>% 
    rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
    gather("grp","lb",-alg,-measure)
  acc_ub <- acc %>% select(alg, measure, contains("ub")) %>%
    rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
    gather("grp","ub",-alg,-measure)
 
  sens_lb <- sens %>% select(alg, measure, contains("lb")) %>% 
    rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
    gather("grp","lb",-alg,-measure)
  sens_ub <- sens %>% select(alg, measure, contains("ub")) %>%
    rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
    gather("grp","ub",-alg,-measure)
 
  spec_lb <- spec %>% select(alg, measure, contains("lb")) %>% 
    rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
    gather("grp","lb",-alg,-measure)
  spec_ub <- spec %>% select(alg, measure, contains("ub")) %>%
    rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
    gather("grp","ub",-alg,-measure)
 
  ppv_lb <- ppv %>% select(alg, measure, contains("lb")) %>% 
    rename(all=lb_all, g1=lb_g1, g2=lb_g2, g3=lb_g3) %>%
    gather("grp","lb",-alg,-measure)
  ppv_ub <- ppv %>% select(alg, measure, contains("ub")) %>%
    rename(all=ub_all, g1=ub_g1, g2=ub_g2, g3=ub_g3) %>%
    gather("grp","ub",-alg,-measure)
  
  }
  acc <- acc %>% select(alg, measure, all, g1, g2, g3) %>% 
    gather("grp","mean",-alg,-measure) %>% left_join(acc_lb) %>% 
    left_join(acc_ub)
  sens <- sens %>% select(alg, measure, all, g1, g2, g3) %>% 
    gather("grp","mean",-alg,-measure) %>% left_join(sens_lb) %>% 
    left_join(sens_ub)
  spec <- spec %>% select(alg, measure, all, g1, g2, g3) %>% 
    gather("grp","mean",-alg,-measure) %>% left_join(spec_lb) %>% 
    left_join(spec_ub)
  ppv <- ppv %>% select(alg, measure, all, g1, g2, g3) %>% 
    gather("grp","mean",-alg,-measure) %>% left_join(ppv_lb) %>% 
    left_join(ppv_ub)
  
  cm_df <- rbind(acc,sens,spec,ppv)
  return(cm_df)
}
bs_wtd <- cmReadFun("bs_wtd") %>% mutate(method="Weighted Bootstrap")
bs_naive <- cmReadFun("bs_naive_pct") %>% mutate(method="Naive Bootstrap")
sp_naive <- cmReadFun("val_naive_in_sample") %>% mutate(method="Naive Standard Practice")

cm_df_data <- rbind(bs_wtd, bs_naive, sp_naive)
cm_df_data$grp <- recode(cm_df_data$grp, "g1"="Class 1","g2"="Class 2","g3"="Class 3","all"="Average")
cm_df_data$alg <- recode(cm_df_data$alg, "enet"="Elastic Net", "gam"="Generalized Additive",
                             "lasso"="Lasso","mn"="Multinomial Logistic",
                             "rf"="Random Forests", "ridge"="Ridge", "xgb"="Gradient Boosting")
cm_df_data <- cm_df_data %>% mutate(method=factor(method, ordered=T, 
                       levels=c("Naive Standard Practice",
                                "Naive Bootstrap",
                                "Weighted Bootstrap"))) 
# suppress error bars -- super tight -- all <0.05
cm_df_data_avg <- cm_df_data %>% filter(grp=="Average")
ggplot(cm_df_data_avg, aes(x=alg, y=mean*100,shape=alg,colour=method)) +
 # geom_errorbar(aes(x=alg, ymin=lb*100, ymax=ub*100), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(0,100))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  #scale_y_continuous(labels = scales::percent)+
  facet_grid(measure ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        panel.background=element_rect(fill="grey92"),
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("Percent")
# Figure 5: Data analysis classification discrimation: Average accuracy, 
  # sensitivity, specificity, and PPV
ggsave(paste0(data_path,"figures/", "data_cm_avg.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

#### DATA SURVIVAL PROBABILTIES : PT ESTIMATES #### 
survPrFun <- function(time){
  surv_bs_wtd <- read.table(paste0(data_path, "bs_wtd_km_", time, ".txt"), sep=",", header=T) 
  surv_bs_naive <- read.table(paste0(data_path,"bs_naive_pct_km_",time,".txt"),sep=",",header=T)
  surv_naive_sample <- read.table(paste0(data_path,"km_naive_pred_in_sample_",time,".txt"),sep=",",header=T)
  
  surv_bs_wtd_lb <- surv_bs_wtd %>% select(alg,contains("lb")) %>% 
    rename(g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class","lb",-alg)
  surv_bs_wtd_ub <- surv_bs_wtd %>% select(alg, contains("ub")) %>%
    rename(g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class","ub",-alg)
  surv_bs_wtd <- surv_bs_wtd %>% select(alg, g1,g2,g3) %>% gather("class","mean",-alg)
  surv_bs_wtd <- surv_bs_wtd %>% left_join(surv_bs_wtd_lb) %>% 
    left_join(surv_bs_wtd_ub) %>% mutate(method="Weighted Bootstrap")
  
  surv_bs_naive_lb <- surv_bs_naive %>% select(alg,contains("lb")) %>% 
    rename(g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class","lb",-alg)
  surv_bs_naive_ub <- surv_bs_naive %>% select(alg, contains("ub")) %>%
    rename(g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class","ub",-alg)
  surv_bs_naive <- surv_bs_naive %>% select(alg, g1,g2,g3) %>% gather("class","mean",-alg)
  surv_bs_naive <- surv_bs_naive %>% left_join(surv_bs_naive_lb) %>% 
    left_join(surv_bs_naive_ub) %>% mutate(method="Naive Bootstrap")
  
  surv_naive_sample_lb <- surv_naive_sample %>% select(alg,contains("lb")) %>% 
    rename(g1=lb_g1, g2=lb_g2, g3=lb_g3) %>% gather("class","lb",-alg)
  surv_naive_sample_ub <- surv_naive_sample %>% select(alg, contains("ub")) %>%
    rename(g1=ub_g1, g2=ub_g2, g3=ub_g3) %>% gather("class","ub",-alg)
  surv_naive_sample <- surv_naive_sample %>% select(alg, g1,g2,g3) %>% gather("class","mean",-alg)
  surv_naive_sample <- surv_naive_sample %>% left_join(surv_naive_sample_lb) %>% 
    left_join(surv_naive_sample_ub) %>% mutate(method="Naive Standard Practice")
  
  # rbind
  surv_plot <- rbind(surv_bs_wtd, surv_bs_naive, surv_naive_sample)
  
  surv_plot$alg <- recode(surv_plot$alg, "enet"="Elastic Net", "gam"="Generalized Additive",
                               "lasso"="Lasso","mn"="Multinomial Logistic",
                               "rf"="Random Forests", "ridge"="Ridge", "xgb"="Gradient Boosting")
  surv_plot$class <- recode(surv_plot$class, "g1"="Stage I/II","g2"="Stage III","g3"="Stage IV")
  surv_plot <- surv_plot %>% mutate(method=factor(method, ordered=T, 
                                                            levels=c("Naive Standard Practice",
                                                                     "Naive Bootstrap",
                                                                     "Weighted Bootstrap")))
  return(surv_plot)
}

# 90 day survival probabilities 
surv_plot_d90 <-survPrFun("d90")
# surpress CIs -- all <.05
ggplot(surv_plot_d90, aes(x=alg, y=mean,shape=alg,colour=method)) +
  #geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  scale_y_continuous()+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("90-Day Survival Probability")
# Figure C7: Data analysis: 90-day survival by predicted stage
ggsave(paste0(data_path,"figures/", "data_surv_d90.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

# 365-day survival probabilities
surv_plot_d365 <-survPrFun("d365")
# supress CIs less than .05 -- all of them
ggplot(surv_plot_d365, aes(x=alg, y=mean,shape=alg,colour=method)) +
  #geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  #scale_y_continuous(labels = scales::percent)+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("365-Day Survival Probability")
# Figure C8: Data analysis: 365-day survival by predicted stage
ggsave(paste0(data_path,"figures/", "data_surv_d365.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)

# Median Survival
surv_plot_med <-survPrFun("median")
# suppress stage IV CIs -- very small on this scale -- 
surv_plot_med <- surv_plot_med %>% mutate(lb=ifelse(ub-lb<31,NA,lb),
                                          ub=ifelse(is.na(lb),NA,ub))
ggplot(surv_plot_med, aes(x=alg, y=mean,shape=alg,colour=method)) +
  geom_errorbar(aes(x=alg, ymin=lb, ymax=ub), width=.15, position="dodge",colour="black") +
  geom_point(size=2) +   
  coord_cartesian(ylim=c(200,600))+
  scale_color_manual(values=okabe[c(5,6,7)]) +
  scale_shape_manual(values=shapePalette)+
  #scale_y_continuous(labels = scales::percent)+
  facet_grid(class ~ method) + 
  theme(panel.border=element_rect(colour="black", fill=NA),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=15),
        strip.text=element_text(size=15),
        panel.background=element_rect(fill="grey92"),
        strip.background=element_rect(colour="black",size=.9,fill=okabe[8]),
        legend.position="none",
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=13, colour="black"),
        axis.text.x=element_text(angle=90,size=13,colour="black",hjust=1,vjust=.25)) +
  ylab("Median Survival (Days)")
# Figure C9: Data analysis: Median survival by predicted stage
ggsave(paste0(data_path,"figures/", "data_surv_med.tiff"), width=8.25, height=9.5, units=c("in"), dpi=500)



