# Check multi-collinearity for liner mixed models ------------------------------
library(car)
library(lmerTest)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("03_LMM_Behaviors_all.RData")
load("04_APC_148.RData")
load("05_Normative.RData")
load("06_longitudinal.RData")

# Behaviors
vars <- c(neurocognition, cbcl, pps, upps, bis)

lmm_vif_behavior <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(vars, df) {
    cat("N =", nrow(df), vars, "\n")
    
    fit <- lmer(df[[vars]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
                  race_ethnicity + family_income + parent_edu + 
                  (1|site_id_l/rel_family_id), data = df)
    
    GVIF <- vif(fit)
    
    return(GVIF)
  }
  
  lmm_behavior <- sapply(vars, lmms, data)
  
  results <- t(lmm_behavior)
  return(results)
}

# Baseline
lmm_behaviors_BL_13 <- lmm_vif_behavior(data_BL, c("Persistently Low", "Increasing"))
lmm_behaviors_BL_12 <- lmm_vif_behavior(data_BL, c("Persistently Low", "Decreasing"))
lmm_behaviors_BL_23 <- lmm_vif_behavior(data_BL, c("Decreasing", "Increasing"))

# 2YFU
lmm_behaviors_2Y_13 <- lmm_vif_behavior(data_2Y, c("Persistently Low", "Increasing"))
lmm_behaviors_2Y_12 <- lmm_vif_behavior(data_2Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_2Y_23 <- lmm_vif_behavior(data_2Y, c("Decreasing", "Increasing"))

# 4YFU
lmm_behaviors_4Y_13 <- lmm_vif_behavior(data_4Y, c("Persistently Low", "Increasing"))
lmm_behaviors_4Y_12 <- lmm_vif_behavior(data_4Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_4Y_23 <- lmm_vif_behavior(data_4Y, c("Decreasing", "Increasing"))

# development change score (Baseline to 2YFU)
lmm_behaviors_BL2Y_13 <- lmm_vif_behavior(data_changes_BL2Y, c("Persistently Low", "Increasing"))
lmm_behaviors_BL2Y_12 <- lmm_vif_behavior(data_changes_BL2Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_BL2Y_23 <- lmm_vif_behavior(data_changes_BL2Y, c("Decreasing", "Increasing"))

# development change score (Baseline to 4YFU)
lmm_behaviors_BL4Y_13 <- lmm_vif_behavior(data_changes_BL4Y, c("Persistently Low", "Increasing"))
lmm_behaviors_BL4Y_12 <- lmm_vif_behavior(data_changes_BL4Y, c("Persistently Low", "Decreasing"))
lmm_behaviors_BL4Y_23 <- lmm_vif_behavior(data_changes_BL4Y, c("Decreasing", "Increasing"))

lmm_behaviors_all_gvif <- rbind(lmm_behaviors_BL_13, lmm_behaviors_BL_12, lmm_behaviors_BL_23,
                                lmm_behaviors_2Y_13, lmm_behaviors_2Y_12, lmm_behaviors_4Y_23,
                                lmm_behaviors_BL2Y_13, lmm_behaviors_BL2Y_12, lmm_behaviors_BL2Y_23,
                                lmm_behaviors_BL4Y_13, lmm_behaviors_BL4Y_12, lmm_behaviors_BL4Y_23) %>%
  as.data.frame()


# Environments -----------------------------------------------------------------
environments <- c(environments[-c(3,4)], peers)
environments[1:8] <- environments[c(1:2, 6:8, 3:5)]
environments[9:17] <- environments[c(9:14, 17, 15, 16)]

lmm_vif_environment <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(vars, df) {
    cat("N =", nrow(df), vars, "\n")
    
    fit <- lmer(df[[vars]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
                  race_ethnicity + family_income + parent_edu + 
                  (1|site_id_l/rel_family_id), data = df)
    
    GVIF <- vif(fit)
    
    return(GVIF)
  }
  
  lmm_environments <- sapply(environments, lmms, data)
  results <- t(lmm_environments)
  
  return(results)
}

# 2YFU
lmm_environments_2Y_13 <- lmm_vif_environment(longitudinal_2Y, c("Persistently Low", "Increasing"))
lmm_environments_2Y_12 <- lmm_vif_environment(longitudinal_2Y, c("Persistently Low", "Decreasing"))
lmm_environments_2Y_23 <- lmm_vif_environment(longitudinal_2Y, c("Decreasing", "Increasing"))

# 4YFU
lmm_environments_4Y_13 <- lmm_vif_environment(longitudinal_4Y, c("Persistently Low", "Increasing"))
lmm_environments_4Y_12 <- lmm_vif_environment(longitudinal_4Y, c("Persistently Low", "Decreasing"))
lmm_environments_4Y_23 <- lmm_vif_environment(longitudinal_4Y, c("Decreasing", "Increasing"))

# development change
lmm_environments_2Y4Y_13 <- lmm_vif_environment(data_changes_2Y4Y, c("Persistently Low", "Increasing"))
lmm_environments_2Y4Y_12 <- lmm_vif_environment(data_changes_2Y4Y, c("Persistently Low", "Decreasing"))
lmm_environments_2Y4Y_23 <- lmm_vif_environment(data_changes_2Y4Y, c("Decreasing", "Increasing"))

lmm_environments_all_gvif <- rbind(
  lmm_environments_2Y_13, lmm_environments_2Y_13, lmm_environments_2Y_23,
  lmm_environments_4Y_13, lmm_environments_4Y_12, lmm_environments_4Y_23,
  lmm_environments_2Y4Y_13, lmm_environments_2Y4Y_12, lmm_environments_2Y4Y_23) %>%
  as.data.frame()


# APC --------------------------------------------------------------------------
smri_volume <- paste0("mrisdp_", 454:601)
smri_area <- paste0("mrisdp_", 303:450)
smri_thick <- paste0("mrisdp_", 1:148)
smri <- c("mrisdp_604", "mrisdp_453", "mrisdp_151")
MRI_data <- c(smri_volume, smri_area, smri_thick, smri_subvolume, smri)

lmm_vif_apc <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(regions, df) {
    cat("N =", nrow(df), regions, "\n")
    
    mod1 <- "df[[regions]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
              race_ethnicity + family_income + parent_edu + site_id_l + 
              smri_vol_scs_intracranialv + 
              (1|mri_info_deviceserialnumber/rel_family_id)"
    
    fit <- lmer(as.formula(mod1), data = df)
    GVIF <- vif(fit)[, 1]
    
    return(GVIF)
  }
  lmm_volume <- sapply(smri_volume, lmms, data)
  lmm_area <- sapply(smri_area, lmms, data)
  lmm_thickness <- sapply(smri_thick, lmms, data)
  lmm_subvolume <- sapply(smri_subvolume, lmms, data)
  lmm_smri <- sapply(smri, lmms, data)
  
  results <- t(cbind(lmm_volume, lmm_area, lmm_thickness, lmm_subvolume, lmm_smri))

  return(results)
}
# BL2Y
lmm_APC_13_BL2Y <- lmm_vif_apc(APC_BL2Y_148, c("Persistently Low", "Increasing"))
lmm_APC_12_BL2Y <- lmm_vif_apc(APC_BL2Y_148, c("Persistently Low", "Decreasing"))
lmm_APC_23_BL2Y <- lmm_vif_apc(APC_BL2Y_148, c("Decreasing", "Increasing"))

# BL4Y
lmm_APC_13_BL4Y <- lmm_vif_apc(APC_BL4Y_148, c("Persistently Low", "Increasing"))
lmm_APC_12_BL4Y <- lmm_vif_apc(APC_BL4Y_148, c("Persistently Low", "Decreasing"))
lmm_APC_23_BL4Y <- lmm_vif_apc(APC_BL4Y_148, c("Decreasing", "Increasing"))

lmm_APC_all_gvif <- rbind(
  lmm_APC_13_BL2Y, lmm_APC_12_BL2Y, lmm_APC_23_BL2Y,
  lmm_APC_13_BL4Y, lmm_APC_12_BL4Y, lmm_APC_23_BL4Y) %>%
  as.data.frame()


# Normative model --------------------------------------------------------------
lmm_vif_normative <- function(data, select_class) {
  data <- filter(data, screen_class %in% select_class)
  
  lmms <- function(regions, df) {
    cat("N =", nrow(df), regions, "\n")
    
    mod1 <- "df[[regions]] ~ screen_class + interview_age + sex + ehi_y_ss_scoreb + 
              race_ethnicity + family_income + parent_edu + site_id_l + 
              smri_vol_scs_intracranialv + 
              (1|mri_info_deviceserialnumber/rel_family_id)"
    
    fit <- lmer(as.formula(mod1), data = df)
    GVIF <- vif(fit)[, 1]
    
    return(GVIF)
  }
  lmm_thickness <- sapply(names(Z_4Y)[3:150], lmms, data)
  lmm_thickness_mean <- sapply(names(Z_4Y)[151:152], lmms, data)

  results <- t(cbind(lmm_thickness, lmm_thickness_mean))
  
  return(results)
}
# Baseline 
lmm_normative_13_BL <- lmm_vif_normative(normative_BL, c("Persistently Low", "Increasing"))
lmm_normative_12_BL <- lmm_vif_normative(normative_BL, c("Persistently Low", "Decreasing"))
lmm_normative_23_BL <- lmm_vif_normative(normative_BL, c("Increasing", "Decreasing"))

# 2YFU
lmm_normative_13_2Y <- lmm_vif_normative(normative_2Y, c("Persistently Low", "Increasing"))
lmm_normative_12_2Y <- lmm_vif_normative(normative_2Y, c("Persistently Low", "Decreasing"))
lmm_normative_23_2Y <- lmm_vif_normative(normative_2Y, c("Increasing", "Decreasing"))

# 4YFU
lmm_normative_13_4Y <- lmm_vif_normative(normative_4Y, c("Persistently Low", "Increasing"))
lmm_normative_12_4Y <- lmm_vif_normative(normative_4Y, c("Persistently Low", "Decreasing"))
lmm_normative_23_4Y <- lmm_vif_normative(normative_4Y, c("Increasing", "Decreasing"))

lmm_environments_all_gvif <- rbind(
  lmm_normative_13_BL, lmm_normative_12_BL, lmm_normative_23_BL,
  lmm_normative_13_2Y, lmm_normative_12_2Y, lmm_normative_23_2Y,
  lmm_normative_13_4Y, lmm_normative_12_4Y, lmm_normative_23_4Y) %>%
  as.data.frame()

# summary results --------------------------------------------------------------
lmm_behaviors_all_gvif$site_id_l <- NA
lmm_behaviors_all_gvif$smri_vol_scs_intracranialv <- NA
lmm_environments_all_gvif$site_id_l <- NA
lmm_environments_all_gvif$smri_vol_scs_intracranialv <- NA

lmm_all_gvif <- rbind(lmm_behaviors_all_gvif, lmm_environments_all_gvif, 
                      lmm_APC_12_BL4Y, lmm_normative_13_BL)

names(lmm_all_gvif) <- c("SMA\nTrajectories", "Age", "Sex", "Handedness", "Race",
                         "Family\nIncome", "Parental\nEducation", 
                         "Study\nSites", "Intracranial\nVolume")

save.image("LMM_GVIF.RData")