# Development changes of environmental factors between baseline and 4YFU
# Bilateral change score models between 2YFU and 4YFU
library(lavaan)
library(lmerTest)
library(tidyverse)

setwd("H:/ABCD/Release5.1/results/25_ABCD_SMA/")
load("02_results_lgca.RData")

# SMA
screens <- c("stq_y_ss_weekday", "screen3_wkdy_y", "screen3b_wkdy_y",
             "screen4_wkdy_y", "screen5_wkdy_y", "screen_wkdy_y")

# environment factors
environments <- environments[-c(3,4)]
environments[1:8] <- environments[c(1:2, 6:8, 3:5)]
environments[9:17] <- environments[c(9:14, 17, 15, 16)]

# select data
longitudinal_2Y <- select(data_2Y, src_subject_id, eventname, screen_class, 
                          all_of(screens), all_of(covarites), 
                          all_of(environments)) %>%
  filter(complete.cases(.))
  
longitudinal_4Y <- select(data_4Y, src_subject_id, eventname, screen_class, 
                          all_of(screens), all_of(covarites), 
                          all_of(environments)) %>%
  filter(complete.cases(.))

common_ID <- list(longitudinal_2Y$src_subject_id, longitudinal_4Y$src_subject_id) %>%
  reduce(intersect)

longitudinal_2Y <- filter(longitudinal_2Y, src_subject_id %in% common_ID)
longitudinal_4Y <- filter(longitudinal_4Y, src_subject_id %in% common_ID)

select(data_2Y, interview_age) %>%
  summarise_all(list(mean, sd))

select(data_4Y, interview_age) %>%
  summarise_all(list(mean, sd))

environment_labels <- c("Family Conflict", "Parental Monitoring",
                        "Neighborhood Crime", "Neighborhood Safety",
                        "Cyberbully", 
                        "School Environment", "School Involvement", "School Disengagement", 
                        "Relational Victimization", "Reputational Aggression",
                        "Reputational Victimization", "Overt Aggression",
                        "Overt Victimization", "Relational Aggression",
                        "Peer Network Health",
                        "Prosocial Peers", "Delinquent Peers")

# Invert variables -------------------------------------------------------------
# Define a function to invert specified variables (higher values indicate worse outcomes)
invert_vars <- function(vars, data) {
  # Check if the variable exists in the data
  if (!vars %in% colnames(data)) {
    stop(paste("Variable", vars, "not found in the data."))
  }
  # Invert the variable: subtract each value from the maximum value in the column
  data[[vars]] <- max(data[[vars]], na.rm = TRUE) - data[[vars]]
  return(data)
}

# List of variables that need to be inverted
vars_need_invert <- c("pmq_y_ss_mean",  "neighborhood_crime_y", 
                      "pbp_ss_prosocial_peers", "srpf_y_ss_ses", 
                      "srpf_y_ss_iiss")

# Apply the function 
longitudinal_2Y <- Reduce(function(data, var) invert_vars(var, data), 
                          vars_need_invert, init = longitudinal_2Y)
longitudinal_4Y <- Reduce(function(data, var) invert_vars(var, data), 
                          vars_need_invert, init = longitudinal_4Y)

# changes ----------------------------------------------------------------------
development_scores <- function(x, data1, data2) {
  changes <- data1[[x]] - data2[[x]]
  changes_score <- lm(changes ~ data2[[x]]) %>%
    residuals()
  return(changes_score)
}

data_changes_2Y4Y <- sapply(environments, development_scores, longitudinal_4Y, longitudinal_2Y) %>%
  cbind(select(longitudinal_2Y, src_subject_id, screen_class, all_of(covarites)))

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
  
  lmm_environments <- sapply(environments, lmms, data)
  
  results <- list(
    "Environments" = lmm_environments
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

# 2YFU
lmm_environments_2Y_13 <- lmm_modality(longitudinal_2Y, c("Persistently Low", "Increasing"))
lmm_environments_2Y_13_fdr <- lapply(lmm_environments_2Y_13, show_fdr)

lmm_environments_2Y_12 <- lmm_modality(longitudinal_2Y, c("Persistently Low", "Decreasing"))
lmm_environments_2Y_12_fdr <- lapply(lmm_environments_2Y_12, show_fdr)

lmm_environments_2Y_23 <- lmm_modality(longitudinal_2Y, c("Decreasing", "Increasing"))
lmm_environments_2Y_23_fdr <- lapply(lmm_environments_2Y_23, show_fdr)

# 4YFU
lmm_environments_4Y_13 <- lmm_modality(longitudinal_4Y, c("Persistently Low", "Increasing"))
lmm_environments_4Y_13_fdr <- lapply(lmm_environments_4Y_13, show_fdr)

lmm_environments_4Y_12 <- lmm_modality(longitudinal_4Y, c("Persistently Low", "Decreasing"))
lmm_environments_4Y_12_fdr <- lapply(lmm_environments_4Y_12, show_fdr)

lmm_environments_4Y_23 <- lmm_modality(longitudinal_4Y, c("Decreasing", "Increasing"))
lmm_environments_4Y_23_fdr <- lapply(lmm_environments_4Y_23, show_fdr)

# development change
lmm_environments_2Y4Y_13 <- lmm_modality(data_changes_2Y4Y, c("Persistently Low", "Increasing"))
lmm_environments_2Y4Y_13_fdr <- lapply(lmm_environments_2Y4Y_13, show_fdr)

lmm_environments_2Y4Y_12 <- lmm_modality(data_changes_2Y4Y, c("Persistently Low", "Decreasing"))
lmm_environments_2Y4Y_12_fdr <- lapply(lmm_environments_2Y4Y_12, show_fdr)

lmm_environments_2Y4Y_23 <- lmm_modality(data_changes_2Y4Y, c("Decreasing", "Increasing"))
lmm_environments_2Y4Y_23_fdr <- lapply(lmm_environments_2Y4Y_23, show_fdr)

# Bivariate change scores models -----------------------------------------------
blcs_func <- function(vars) {
  blcs <- function(x, vars) {
    
    # regressed out covariates
    residual_value <- function(y, data) {
      fit <- lmer(data[[y]] ~ interview_age + sex + race_ethnicity + 
                    ehi_y_ss_scoreb +  family_income  + parent_edu + 
                    (1|site_id_l/rel_family_id), 
                  data = data)
      return(residuals(fit))
    }
    
    residual_2Y_stqweekday <- residual_value(vars, longitudinal_2Y)
    residual_4Y_stqweekday <- residual_value(vars, longitudinal_4Y)
    
    residual_2Y_x <- residual_value(x, longitudinal_2Y)
    residual_4Y_x <- residual_value(x, longitudinal_4Y)
    
    residualed_2Y <- data.frame(
      "src_subject_id" = longitudinal_2Y$src_subject_id,
      "stqweekday_2Y" = residual_2Y_stqweekday,
      "x_2Y" = residual_2Y_x
    )
    residualed_4Y <- data.frame(
      "src_subject_id" = longitudinal_4Y$src_subject_id,
      "stqweekday_4Y" = residual_4Y_stqweekday,
      "x_4Y" = residual_4Y_x
    )
    
    blcs_data <- inner_join(residualed_2Y, residualed_4Y, by = "src_subject_id") %>%
      select(-src_subject_id)
    
    # the Bivariate Latent Change Score model --------------------------
    BLCS <- '
    stqweekday_4Y ~ 1*stqweekday_2Y
    dstqweekday =~ 1*stqweekday_4Y
    dstqweekday ~ 1
    stqweekday_2Y ~ 1
    stqweekday_4Y ~ 0*1

    x_4Y ~ 1*x_2Y
    dx =~ 1*x_4Y
    x_4Y ~ 0*1
    x_4Y ~~ 0*x_4Y

    dstqweekday ~~ dstqweekday
    stqweekday_2Y ~~ stqweekday_2Y
    stqweekday_4Y ~~ 0*stqweekday_4Y

    dx ~ 1
    x_2Y ~ 1
    dx ~~ dx
    x_2Y ~~ x_2Y

    dx ~ stqweekday_2Y + x_2Y
    dstqweekday ~ x_2Y + stqweekday_2Y

    stqweekday_2Y ~~ x_2Y
    dstqweekday ~~ dx
  '
    
    fits <- lavaan(BLCS, data = blcs_data, estimator = "mlr", fixed.x = FALSE)
    result <- standardizedsolution(fits, level = 0.95)[17:20, ]
    return(result)
  }
  
  blcs_all <- lapply(environments, blcs, vars)
  names(blcs_all) <- environment_labels
  
  
  blcs_results <- function(index, data) {
    df <- data.frame()
    
    for(i in seq(data)) {
      df <- rbind(df, unlist(data[[i]][index, ]))
    }
    colnames(df) <- colnames(blcs_all[[1]])
    
    # FDR corrections
    df$fdr <- p.adjust(df$pvalue, method = "fdr")
    
    df <- select(df, est.std, ci.lower, ci.upper, se, pvalue, fdr) %>%
      summarise_all(as.numeric) %>%
      mutate(environment = environment_labels)
    
    df[, 1:4] <- apply(df[, 1:4], 2, round, 3)
    df$screens <- vars
    # df <- filter(df, fdr < 0.05) %>%
    #   arrange(est.std)
    
    return(df)
  }
  BLstqweekday_4Yenvironments <- blcs_results(1, blcs_all)
  BLenvironments_4Ystqweekday <- blcs_results(3, blcs_all)
  
  cat(vars)
  
  return(list(
    "A" = BLstqweekday_4Yenvironments, 
    "B" = BLenvironments_4Ystqweekday))
  
}
 
blcs_results_1 <- blcs_func(screens[1])  # overall SMA
blcs_results_2 <- blcs_func(screens[2])  # playing video games (single-player)
blcs_results_3 <- blcs_func(screens[3])  # playing video games (multi-player)
blcs_results_4 <- blcs_func(screens[4])  # texting on digital devices
blcs_results_5 <- blcs_func(screens[5])  # visiting social network sites
blcs_results_6 <- blcs_func(screens[6])  # video chatting (no significant results)

save.image("06_longitudinal.RData")