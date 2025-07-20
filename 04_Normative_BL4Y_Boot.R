# Bootstrap approach for balancing SMA trajectory sample sizes
# Normative Models at 4YFU
library(lmerTest)
library(doParallel)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("05_Normative_BL4Y.RData")

# boot -------------------------------------------------------------------------
table(normative_4Y$screen_class)

# Extract the average participants
Low_4Y <- which(normative_4Y$screen_class == "Persistently Low")
n_4Y <- round((sum(normative_4Y$screen_class == "Increasing") + sum(normative_4Y$screen_class == "Decreasing")) / 2)

# 1000 times
nboot <- 1000

set.seed(123)
Low_4Y_Boot <- replicate(nboot, sample(Low_4Y, n_4Y, replace = FALSE))

# LMM boot ---------------------------------------------------------------------
lmm_modality_boot <- function(boot_index, boot_mat, data, select_class) {
  Increasing <- which(data$screen_class == "Increasing")
  Decreasing <- which(data$screen_class == "Decreasing")
  Low <- boot_mat[, boot_index]
  
  include_index <- c(Low, Increasing, Decreasing)
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

# set do parallel --------------------------------------------------------------
cl <- makeCluster(14, type = 'PSOCK')
clusterExport(cl, c("normative_4Y", "Low_4Y_Boot"))
registerDoParallel(cl)

# Persistently low vs. Decreasing
lmm_normative_12_4Y_boot <- foreach(boot_seq = seq(nboot), 
                                    .verbose = TRUE,
                                    .packages = c("tidyverse", "lmerTest")) %dopar% 
  lmm_modality_boot(boot_seq, Low_4Y_Boot, normative_4Y, 
                    c("Persistently Low", "Decreasing"))

# Persistently low vs. Increasing
lmm_normative_13_4Y_boot <- foreach(boot_seq = seq(nboot),
                                    .verbose = TRUE,
                                    .packages = c("tidyverse", "lmerTest")) %dopar% 
  lmm_modality_boot(boot_seq, Low_4Y_Boot, normative_4Y, 
                    c("Persistently Low", "Increasing"))

save.image("05_Normative_4Y_Boot.RData")