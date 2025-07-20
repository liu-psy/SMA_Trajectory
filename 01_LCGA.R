# The latent class linear mixed models
# SMA from baseline to 4YFU
library(lcmm)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA")
load("01_sample_selection.RData")

# LCGA -------------------------------------------------------------------------
vars <- c("src_subject_id", "eventname", "stq_y_ss_weekday", "interview_age",
          "sex")

screen_BL <- select(data_BL, all_of(vars))
screen_1Y <- select(data_1Y, all_of(vars))
screen_2Y <- select(data_2Y, all_of(vars))
screen_3Y <- select(data_3Y, all_of(vars))
screen_4Y <- select(data_4Y, all_of(vars))

analysis_data <- rbind(screen_BL, screen_1Y, screen_2Y, screen_3Y, screen_4Y) %>%
  mutate("ID" = as.numeric(as.factor(src_subject_id)))

# gridsearch -------------------------------------------------------------------
lcga1 <- hlme(stq_y_ss_weekday ~ interview_age + sex,
              random = ~ interview_age + sex, 
              subject = "ID",
              ng = 1, 
              data = analysis_data)

# class 2-8
gridsearch_func <- function(nclass) {
  results <- gridsearch(rep = 500, maxiter = 50, minit = lcga1, cl = 16, 
                        m = hlme(stq_y_ss_weekday ~ interview_age + sex,
                                 random = ~ interview_age + sex,
                                 mixture = ~ interview_age + sex,
                                 subject = "ID",
                                 ng = nclass, 
                                 data = analysis_data))
  return(results)
}
lcga2 <- gridsearch_func(nclass = 2)
lcga3 <- gridsearch_func(nclass = 3)
lcga4 <- gridsearch_func(nclass = 4)
lcga5 <- gridsearch_func(nclass = 5)
lcga6 <- gridsearch_func(nclass = 6)
lcga7 <- gridsearch_func(nclass = 7)
lcga8 <- gridsearch_func(nclass = 8)

summarytable(lcga1, lcga2, lcga3, lcga4, lcga5, lcga6, lcga7, lcga8,
             which = c("G", "loglik", "conv", "npm", "AIC", "BIC", "SABIC", 
                       "entropy","ICL", "%class"))

# 3 trajectories were selected
screen_class <- factor(lcga3$pprob$class, levels = 1:3,  
                       labels = c("Persistently Low", "Decreasing", "Increasing"))

# Statistics tests -------------------------------------------------------------
data_BL$screen_class <- screen_class

# Kruskal-Wallis test
kruskal.test(interview_age ~ screen_class, data = data_BL)
kruskal.test(family_income ~ screen_class, data = data_BL)
kruskal.test(parent_edu ~ screen_class, data = data_BL)
# Chi-squared test
chisq.test(table_sex)
chisq.test(table_race)
chisq.test(table_handedness)