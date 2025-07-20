# Normative models at 4YFU
library(lmerTest)
library(neuroCombat)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("04_APC_BL4Y.RData")

# load results from normative model 
Z_BL <- read_csv("H:/ABCD/Release5.1/results/normative_elife/normative_BL_148/Z.csv")
Z_4Y <- read_csv("H:/ABCD/Release5.1/results/normative_elife/normative_4Y_148/Z.csv")

normative_BL <- select(data_BL, src_subject_id, screen_class, eventname, 
                       all_of(covarites), all_of(mri_qc), 
                       smri_vol_scs_intracranialv) %>%
  inner_join(Z_BL, by = c("src_subject_id", "eventname")) %>%
  filter(imgincl_t1w_include == 1) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

normative_4Y <- select(data_4Y, src_subject_id, screen_class, eventname, 
                       all_of(covarites), all_of(mri_qc), 
                       smri_vol_scs_intracranialv) %>%
  inner_join(Z_4Y, by = c("src_subject_id", "eventname")) %>%
  filter(imgincl_t1w_include == 1) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

# LMM --------------------------------------------------------------------------
lmm_modality <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(regions, df) {
    cat("N =", nrow(df), regions, "\n")
    
    mod1 <- "df[[regions]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
              race_ethnicity + family_income + parent_edu + site_id_l + 
              smri_vol_scs_intracranialv + 
              (1|mri_info_deviceserialnumber/rel_family_id)"
    
    fit <- lmer(as.formula(mod1), data = df)
    
    beta_CI <- round(confint(fit, method = "Wald")[5, ], 3)
    beta <- round(coef(summary(fit))[2, 1], 2)
    beta_se <- round(coef(summary(fit))[2, 2], 3)
    t <- round(coef(summary(fit))[2, 4], 2)
    p <- coef(summary(fit))[2, 5]
    
    results <- c(t, beta, beta_se, beta_CI, p)
    names(results) <- c("t", "beta", "se", "2.5%", "97.5%", "p")
    
    return(results)
  }
  lmm_thickness <- sapply(names(Z_4Y)[3:150], lmms, data)
  lmm_thickness_mean <- sapply(names(Z_4Y)[151:152], lmms, data)
  
  results <- list(
    "Thickness" = lmm_thickness,
    "Thickness_mean" = lmm_thickness_mean
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
  lmm_results <- as_tibble(lmm_results, rownames = NA)
  lmm_results <- lmm_results[, as.vector(lmm_results["p_fdr", ] < 0.05)]
  return(lmm_results)
}

# Baseline 
lmm_normative_13_BL <- lmm_modality(normative_BL, c("Persistently Low", "Increasing"))
lmm_normative_13_BL_fdr <- lapply(lmm_normative_13_BL, show_fdr)

lmm_normative_12_BL <- lmm_modality(normative_BL, c("Persistently Low", "Decreasing"))
lmm_normative_12_BL_fdr <- lapply(lmm_normative_12_BL, show_fdr)

lmm_normative_23_BL <- lmm_modality(normative_BL, c("Increasing", "Decreasing"))
lmm_normative_23_BL_fdr <- lapply(lmm_normative_23_BL, show_fdr)

# 4YFU
lmm_normative_13_4Y <- lmm_modality(normative_4Y, c("Persistently Low", "Increasing"))
lmm_normative_13_4Y_fdr <- lapply(lmm_normative_13_4Y, show_fdr)

lmm_normative_12_4Y <- lmm_modality(normative_4Y, c("Persistently Low", "Decreasing"))
lmm_normative_12_4Y_fdr <- lapply(lmm_normative_12_4Y, show_fdr)

lmm_normative_23_4Y <- lmm_modality(normative_4Y, c("Increasing", "Decreasing"))
lmm_normative_23_4Y_fdr <- lapply(lmm_normative_23_4Y, show_fdr)

# output t-maps
df_tmaps <- data.frame(
  # Normative
  "Normative_12" = lmm_normative_12_4Y$Thickness[1, ],
  "Normative_13" = lmm_normative_13_4Y$Thickness[1, ],
  # APC
  "APC_12" = lmm_APC_12_BL4Y$Thickness[1, ],
  "APC_13" = lmm_APC_13_BL4Y$Thickness[1, ]
  )
write_csv(df_tmaps, "APC_normative_tmaps.csv")

# save results
save.image("05_Normative_BL4Y.RData")