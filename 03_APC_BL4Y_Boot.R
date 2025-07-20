# Bootstrap approach for balancing SMA trajectory sample sizes
# Brain annualized percent change between baseline and 4YFU
library(lmerTest)
library(doParallel)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("04_APC_BL4Y.RData")

# boot -------------------------------------------------------------------------
table(apc_BL4Y_BL$screen_class)

# Extract the average participants
Low_BL4Y <- which(APC_BL4Y_148$screen_class == "Persistently Low")
BL4Y_n <- round((sum(APC_BL4Y_148$screen_class == "Increasing") + sum(APC_BL4Y_148$screen_class == "Decreasing")) / 2)

# 1000 times
nboot <- 1000

set.seed(123)
Low_BL4Y_Boot <- replicate(nboot, sample(Low_BL4Y, BL4Y_n, replace = FALSE))

# LMM boot ---------------------------------------------------------------------
lmm_modality_boot <- function(boot_index, boot_mat, data, select_class) {
  Increasing_BL2Y <- which(data$screen_class == "Increasing")
  Decreasing_BL2Y <- which(data$screen_class == "Decreasing")
  Low_BL2Y <- boot_mat[, boot_index]
  
  include_index <- c(Low_BL2Y, Increasing_BL2Y, Decreasing_BL2Y)
  data <- data[include_index, ]
  
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
  
  lmm_thickness <- sapply(smri_thick, lmms, data)
  lmm_smri <- sapply(smri, lmms, data)
  
  results <- list(
    "Thickness" = lmm_thickness,
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

# set do parallel --------------------------------------------------------------
cl <- makeCluster(14, type = 'PSOCK')
clusterExport(cl, c("APC_BL4Y_148",  "Low_BL4Y_Boot"))
registerDoParallel(cl)

# BL4Y
lmm_APC_12_BL4Y_boot <- foreach(boot_seq = seq(nboot), 
                                .verbose = TRUE,
                                .packages = c("tidyverse", "lmerTest")) %dopar% 
  lmm_modality_boot(boot_seq, Low_BL4Y_Boot, APC_BL4Y_148, 
                    c("Persistently Low", "Decreasing"))

lmm_APC_13_BL4Y_boot <- foreach(boot_seq = seq(nboot),
                                .verbose = TRUE,
                                .packages = c("tidyverse", "lmerTest")) %dopar% 
  lmm_modality_boot(boot_seq, Low_BL4Y_Boot, APC_BL4Y_148, 
                    c("Persistently Low", "Increasing"))

save.image("04_APC_BL4Y_Boot.RData")