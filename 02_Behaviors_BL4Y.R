# Development change score of behavioral phenotype between baseline and 4YFU
library(lmerTest)
library(openxlsx)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("02_results_lgca.RData")

vars <- c(neurocognition, cbcl, pps, upps, bis)

data_BL <- select(data_BL, src_subject_id, eventname, screen_class, 
                  all_of(covarites), all_of(vars)) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

data_4Y <- select(data_4Y, src_subject_id, eventname, screen_class, 
                  all_of(covarites), all_of(vars)) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

# common subjects
common_ID_BL4Y <- list(data_BL$src_subject_id, data_4Y$src_subject_id) %>%
  reduce(intersect)

data_BL4Y_BL <- filter(data_BL, src_subject_id %in% common_ID_BL4Y)
data_BL4Y_4Y <- filter(data_4Y, src_subject_id %in% common_ID_BL4Y)

select(data_BL, interview_age) %>%
  summarise_all(list(mean, sd))
select(data_4Y, interview_age) %>%
  summarise_all(list(mean, sd))

# development changes ----------------------------------------------------------
development_scores <- function(x, data1, data2) {
  changes <- data1[[x]] - data2[[x]]
  changes_score <- lm(changes ~ data2[[x]]) %>%
    residuals()
  return(changes_score)
}

data_changes_BL4Y <- sapply(vars, development_scores, data_BL4Y_4Y, data_BL4Y_BL) %>%
  cbind(select(data_BL4Y_BL, src_subject_id, screen_class, all_of(covarites)))

# LMM --------------------------------------------------------------------------
lmm_modality <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(vars, df) {
    cat("N =", nrow(df), vars, "\n")
    
    fit <- lmer(df[[vars]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
                  race_ethnicity + family_income + parent_edu + 
                  (1|site_id_l/rel_family_id), data = df)
    
    beta_CI <- round(confint(fit, method = "Wald")[5, ], 3)
    beta <- round(coef(summary(fit))[2, 1], 2)
    beta_se <- round(coef(summary(fit))[2, 2], 3)
    t <- round(coef(summary(fit))[2, 4], 2)
    p <- coef(summary(fit))[2, 5]
    
    results <- c(t, beta, beta_se, beta_CI, p)
    names(results) <- c("t", "beta", "se", "2.5%", "97.5%", "p")
    
    return(results)
  }
  
  lmm_behavior <- sapply(vars, lmms, data)
  
  results <- list(
    "Bahviors" = lmm_behavior
  )
  
  # FDR corrections
  p_correct <- function(lmm_results) {
    p_value <- lmm_results[6, ]
    p_fdr <- p.adjust(p_value, "fdr")
    lmm_results <- rbind(lmm_results, p_fdr)
    return(lmm_results)
  }
  
  results <- lapply(results, p_correct)
  
  return(results)
}

# FDR corrections
show_fdr <- function(lmm_results) {
  lmm_results <- lmm_results[, lmm_results["p_fdr", ] < 0.05]
  return(lmm_results)
}

# Baseline
lmm_behaviors_BL_13 <- lmm_modality(data_BL, c("Persistently Low", "Increasing"))
lmm_behaviors_BL_13_fdr <- lapply(lmm_behaviors_BL_13, show_fdr)

lmm_behaviors_BL_12 <- lmm_modality(data_BL, c("Persistently Low", "Decreasing"))
lmm_behaviors_BL_12_fdr <- lapply(lmm_behaviors_BL_12, show_fdr)

lmm_behaviors_BL_23 <- lmm_modality(data_BL, c("Decreasing", "Increasing"))
lmm_behaviors_BL_23_fdr <- lapply(lmm_behaviors_BL_23, show_fdr)

# 4YFU
lmm_behaviors_4Y_13 <- lmm_modality(data_4Y, c("Persistently Low", "Increasing"))
lmm_behaviors_4Y_13_fdr <- lapply(lmm_behaviors_4Y_13, show_fdr)

lmm_behaviors_4Y_12 <- lmm_modality(data_4Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_4Y_12_fdr <- lapply(lmm_behaviors_4Y_12, show_fdr)

lmm_behaviors_4Y_23 <- lmm_modality(data_4Y, c("Decreasing", "Increasing"))
lmm_behaviors_4Y_23_fdr <- lapply(lmm_behaviors_4Y_23, show_fdr)

# development change score (Baseline to 4YFU)
lmm_behaviors_BL4Y_13 <- lmm_modality(data_changes_BL4Y, c("Persistently Low", "Increasing"))
lmm_behaviors_BL4Y_13_fdr <- lapply(lmm_behaviors_BL4Y_13, show_fdr)

lmm_behaviors_BL4Y_12 <- lmm_modality(data_changes_BL4Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_BL4Y_12_fdr <- lapply(lmm_behaviors_BL4Y_12, show_fdr)

lmm_behaviors_BL4Y_23 <- lmm_modality(data_changes_BL4Y, c("Decreasing", "Increasing"))
lmm_behaviors_BL4Y_23_fdr <- lapply(lmm_behaviors_BL4Y_23, show_fdr)

# save results
save.image("03_LMM_Behaviors_BL4Y.RData")