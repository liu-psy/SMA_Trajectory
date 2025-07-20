# Brain annualized percent change between baseline and 4YFU
library(lmerTest)
library(mediation)
library(neuroCombat)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("02_results_lgca.RData")

mri_destrieux <- read_csv("H:/ABCD/Release5.1/mri_destrieux.csv")
data_BL <- left_join(data_BL, mri_destrieux, by = c("src_subject_id", "eventname"))
data_4Y <- left_join(data_4Y, mri_destrieux, by = c("src_subject_id", "eventname"))

# APCs -------------------------------------------------------------------------
# MRI (destrieux)
smri_volume <- paste0("mrisdp_", 454:601)
smri_area <- paste0("mrisdp_", 303:450)
smri_thick <- paste0("mrisdp_", 1:148)
smri <- c("mrisdp_604", "mrisdp_453", "mrisdp_151")

MRI_data <- c(smri_volume, smri_area, smri_thick, smri_subvolume, smri)
mri_qc <- c("mri_info_deviceserialnumber", "imgincl_t1w_include")

# APCs -------------------------------------------------------------------------
apc_BL <- select(data_BL, src_subject_id, screen_class, eventname, interview_age, 
                 all_of(MRI_data), all_of(mri_qc)) %>%
  filter(imgincl_t1w_include == 1) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

apc_4Y <- select(data_4Y, src_subject_id, screen_class, eventname, interview_age, 
                 all_of(MRI_data), all_of(mri_qc)) %>%
  filter(imgincl_t1w_include == 1) %>%
  filter(complete.cases(.)) %>%
  filter(!duplicated(.))

# common subjects
common_ID_BL4Y <- list(apc_BL$src_subject_id, apc_4Y$src_subject_id) %>%
  reduce(intersect)

apc_BL4Y_BL <- filter(apc_BL, src_subject_id %in% common_ID_BL4Y)
apc_BL4Y_4Y <- filter(apc_4Y, src_subject_id %in% common_ID_BL4Y)

covariate_BL4Y <- select(data_BL, src_subject_id, screen_class, all_of(covarites), 
                         mri_info_deviceserialnumber, smri_vol_scs_intracranialv) %>%
  filter(src_subject_id %in% common_ID_BL4Y)

# 1. Combat --------------------------------------------------------------------
harmonized_data <- function(data, vars) {
  # neurocombat
  neurocombat_data <- neuroCombat(t(select(data, all_of(vars))),
                                  batch = data$mri_info_deviceserialnumber)
  
  neurocombat_data <- t(neurocombat_data$dat.combat) %>%
    as.data.frame()
  
  neurocombat_data <- cbind(neurocombat_data, 
                            select(data, src_subject_id, interview_age))
  
  return(neurocombat_data)
}
harmonized_BL4Y_BL <- harmonized_data(apc_BL4Y_BL, MRI_data)
harmonized_BL4Y_4Y <- harmonized_data(apc_BL4Y_4Y, MRI_data)

# Annualized percent change APCs -----------------------------------------------
APC <- function(index, data1, data2) {
  delta_age <- (data1$interview_age - data2$interview_age)
  a1 <- data1[, index] - data2[, index]
  a2 <- data1[, index] + data2[, index]
  apc <- ((a1/a2) * 100) / delta_age
  
  return(apc)
}
APC_BL4Y_148 <- sapply(MRI_data, APC, harmonized_BL4Y_4Y, harmonized_BL4Y_BL) %>%
  as.data.frame() 
APC_BL4Y_148 <- cbind(covariate_BL4Y, APC_BL4Y_148) %>%
  filter(complete.cases(.))

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
  lmm_volume <- sapply(smri_volume, lmms, data)
  lmm_area <- sapply(smri_area, lmms, data)
  lmm_thickness <- sapply(smri_thick, lmms, data)
  lmm_subvolume <- sapply(smri_subvolume, lmms, data)
  lmm_smri <- sapply(smri, lmms, data)
  
  results <- list(
    "Volume" = lmm_volume,
    "Areas" = lmm_area,
    "Thickness" = lmm_thickness,
    "Subvolume" = lmm_subvolume,
    "Global_sMRI" = lmm_smri
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

# BL4Y
lmm_APC_13_BL4Y <- lmm_modality(APC_BL4Y_148, c("Persistently Low", "Increasing"))
lmm_APC_13_BL4Y_fdr <- lapply(lmm_APC_13_BL4Y, show_fdr)

lmm_APC_12_BL4Y <- lmm_modality(APC_BL4Y_148, c("Persistently Low", "Decreasing"))
lmm_APC_12_BL4Y_fdr <- lapply(lmm_APC_12_BL4Y, show_fdr)

lmm_APC_23_BL4Y <- lmm_modality(APC_BL4Y_148, c("Decreasing", "Increasing"))
lmm_APC_23_BL4Y_fdr <- lapply(lmm_APC_23_BL4Y, show_fdr)

# Mediation analysis -----------------------------------------------------------
mediation_analysis <- function(sma_class, vars) {
  vars_BL <- dplyr::select(data_BL, src_subject_id, vars) %>%
    filter(complete.cases(.))
  vars_4Y <- dplyr::select(data_4Y, src_subject_id, vars) %>%
    filter(complete.cases(.))
  
  vars_change <- inner_join(vars_BL, vars_4Y, by = "src_subject_id")
  
  vars_development <- vars_change[[3]] - vars_change[[2]]
  development_score <- lm(vars_development ~ vars_change[[2]]) %>%
    residuals()
  
  vars_change <- mutate(vars_change, "development_score" = development_score) %>%
    select(src_subject_id, development_score)
  
  
  mediation_data <- inner_join(APC_BL4Y_148, vars_change, by = "src_subject_id") %>%
    filter(screen_class %in% c("Persistently Low", sma_class)) %>%
    mutate(
      "screen_class" = as.numeric(factor(screen_class)),
      "site_id_l" = as.numeric(factor(site_id_l)),
      "mri_info_deviceserialnumber" = as.numeric(factor(mri_info_deviceserialnumber)),
    ) %>%
    dplyr::select(screen_class, development_score, mrisdp_604, mrisdp_151,
                  all_of(covarites), mri_info_deviceserialnumber, 
                  smri_vol_scs_intracranialv) %>%
    filter(complete.cases(.))
  
  print(table(mediation_data$screen_class))
  
  ## mediation analysis
  # mediation model
  mediate_model <- lm(mrisdp_151 ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
                        race_ethnicity + family_income + parent_edu + site_id_l + 
                        mri_info_deviceserialnumber + smri_vol_scs_intracranialv, data = mediation_data)
  # full model
  full_model <- lm(development_score ~ screen_class + mrisdp_151 + interview_age + sex + ehi_y_ss_scoreb + 
                     race_ethnicity + family_income + parent_edu + site_id_l + 
                     mri_info_deviceserialnumber + smri_vol_scs_intracranialv, data = mediation_data)
  # mediation
  mediation_fit <- mediate(mediate_model, full_model, treat = "screen_class", mediator = "mrisdp_151", 
                           boot = TRUE, sims = 10000)
  
  return(mediation_fit)
}

mediation_13_rulebreak <- mediation_analysis("Increasing", "cbcl_scr_syn_rulebreak_r")
summary(mediation_13_rulebreak)

mediation_13_planning <- mediation_analysis("Increasing", "upps_y_ss_lack_of_planning")
summary(mediation_13_planning)

# save results
save.image("04_APC_BL4Y.RData")

# Annualized percent change APCs -----------------------------------------------
# Regressed out baseline value (controlling for baseline brain metrics)
APC <- function(index, data1, data2) {
  delta_age <- (data1$interview_age - data2$interview_age)
  a1 <- data1[, index] - data2[, index]
  a2 <- data1[, index] + data2[, index]
  apc <- ((a1/a2) * 100) / delta_age
  
  # Regressed out baseline value
  apc <- lm(apc ~ data2[, index]) %>%
    residuals()
  
  return(apc)
}