#!/usr/bin/env Rscript
# ============================================================================
# Mixed-Effects Models for RS-Derived Diversity & Taxonomic Diversity
# ============================================================================
# Fits lme4 models linking remote-sensing-derived diversity to taxonomic
# diversity and productivity, following the NASA NPP proposal framework.
#
# Models:
#   1a) Alpha diversity ~ spectral * structural + heterogeneity + (1|domain/site)
#   1b) Beta diversity  ~ spectral * structural + heterogeneity + (1|domain/site)
#   2)  Productivity    ~ functional_div * heterogeneity + (1|domain/site)
#
# Usage:
#   Rscript analysis/run_mixed_models.R
# ============================================================================

library(lme4)
library(lmerTest)  # p-values for lmer
library(MuMIn)     # r.squaredGLMM
library(ggplot2)
library(dplyr)
library(car)       # vif

# ── Paths ──
base_dir <- "E:/neon_lidar/model_results"
fig_dir  <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──
cat("Loading assembled data...\n")
df <- read.csv(file.path(base_dir, "assembled_data.csv"), stringsAsFactors = FALSE)
cat(sprintf("  %d rows, %d columns, %d sites\n", nrow(df), ncol(df), length(unique(df$siteID))))

# Convert grouping variables to factors
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)
df$year   <- as.numeric(df$year)

# ── Standardize predictors (z-score) ──
scale_if_exists <- function(df, col) {
  if (col %in% names(df) && sum(!is.na(df[[col]])) > 5) {
    df[[paste0(col, "_z")]] <- scale(df[[col]])[, 1]
  }
  df
}

pred_cols <- c("rao_q", "spectral_cv", "spectral_shannon",
               "spectral_FRic", "spectral_FDiv",
               "chm_cv", "chm_sd", "Sa_500m", "Sq_500m",
               "func_FRic", "func_FDiv", "func_RaoQ",
               "nlcd_shdi", "cumulative_mean", "trend_mean")

for (col in pred_cols) {
  df <- scale_if_exists(df, col)
}

# ── Helper function ──
fit_and_report <- function(formula_str, data, model_name) {
  cat(sprintf("\n══════════════════════════════════════\n"))
  cat(sprintf("Model: %s\n", model_name))
  cat(sprintf("Formula: %s\n", formula_str))

  tryCatch({
    m <- lmer(as.formula(formula_str), data = data, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))

    cat("\n--- Summary ---\n")
    print(summary(m))

    # R-squared (marginal and conditional)
    r2 <- r.squaredGLMM(m)
    cat(sprintf("\nR2 marginal (fixed):      %.4f\n", r2[1, "R2m"]))
    cat(sprintf("R2 conditional (fixed+RE): %.4f\n", r2[1, "R2c"]))

    # VIF for fixed effects
    cat("\n--- VIF ---\n")
    tryCatch({
      v <- vif(m)
      print(v)
    }, error = function(e) cat("  VIF computation failed (likely interaction term)\n"))

    # Save coefficients
    coef_df <- as.data.frame(coef(summary(m)))
    coef_df$term <- rownames(coef_df)
    coef_df$model <- model_name
    coef_df$R2m <- r2[1, "R2m"]
    coef_df$R2c <- r2[1, "R2c"]

    # Diagnostic plot
    png(file.path(fig_dir, paste0(gsub(" ", "_", model_name), "_diagnostics.png")),
        width = 800, height = 600)
    par(mfrow = c(2, 2))
    plot(fitted(m), resid(m), main = "Residuals vs Fitted",
         xlab = "Fitted", ylab = "Residuals", pch = 16, col = rgb(0,0,0,0.3))
    abline(h = 0, col = "red")
    qqnorm(resid(m), main = "Q-Q Plot", pch = 16, col = rgb(0,0,0,0.3))
    qqline(resid(m), col = "red")
    hist(resid(m), breaks = 30, main = "Residual Distribution", col = "lightblue")
    plot(fitted(m), sqrt(abs(resid(m))), main = "Scale-Location",
         xlab = "Fitted", ylab = "sqrt(|Residuals|)", pch = 16, col = rgb(0,0,0,0.3))
    dev.off()
    cat(sprintf("  Saved: %s_diagnostics.png\n", gsub(" ", "_", model_name)))

    return(coef_df)

  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
}


# ════════════════════════════════════════════════════════════════════════════
# Model 1a: Alpha taxonomic diversity ~ RS diversity
# ════════════════════════════════════════════════════════════════════════════

all_coefs <- list()

# Check which predictors are available
has_spectral <- "rao_q_z" %in% names(df) && sum(!is.na(df$rao_q_z)) > 10
has_structural <- "chm_cv_z" %in% names(df) && sum(!is.na(df$chm_cv_z)) > 10
has_heterogeneity <- "Sa_500m_z" %in% names(df) && sum(!is.na(df$Sa_500m_z)) > 10

if ("alpha_shannon_mean" %in% names(df) && has_spectral && has_structural) {
  # Full interaction model
  f1a <- "alpha_shannon_mean ~ rao_q_z * chm_cv_z + (1|domain/siteID)"
  c1a <- fit_and_report(f1a, df, "M1a Alpha Shannon")
  if (!is.null(c1a)) all_coefs[[length(all_coefs) + 1]] <- c1a

  if (has_heterogeneity) {
    f1a2 <- "alpha_shannon_mean ~ rao_q_z * chm_cv_z + Sa_500m_z + (1|domain/siteID)"
    c1a2 <- fit_and_report(f1a2, df, "M1a Alpha Shannon + Heterogeneity")
    if (!is.null(c1a2)) all_coefs[[length(all_coefs) + 1]] <- c1a2
  }

  # Species richness model
  if ("alpha_richness_mean" %in% names(df)) {
    f1a_r <- "alpha_richness_mean ~ rao_q_z * chm_cv_z + (1|domain/siteID)"
    c1a_r <- fit_and_report(f1a_r, df, "M1a Alpha Richness")
    if (!is.null(c1a_r)) all_coefs[[length(all_coefs) + 1]] <- c1a_r
  }
} else {
  cat("\nSkipping Model 1a: insufficient predictors\n")
  cat(sprintf("  has_spectral=%s, has_structural=%s\n", has_spectral, has_structural))
}


# ════════════════════════════════════════════════════════════════════════════
# Model 1b: Beta taxonomic diversity ~ RS diversity
# ════════════════════════════════════════════════════════════════════════════

if ("bray_mean" %in% names(df) && has_spectral && has_structural) {
  f1b <- "bray_mean ~ rao_q_z * chm_cv_z + (1|domain/siteID)"
  c1b <- fit_and_report(f1b, df, "M1b Beta Bray-Curtis")
  if (!is.null(c1b)) all_coefs[[length(all_coefs) + 1]] <- c1b

  if (has_heterogeneity) {
    f1b2 <- "bray_mean ~ rao_q_z * chm_cv_z + Sa_500m_z + (1|domain/siteID)"
    c1b2 <- fit_and_report(f1b2, df, "M1b Beta Bray-Curtis + Heterogeneity")
    if (!is.null(c1b2)) all_coefs[[length(all_coefs) + 1]] <- c1b2
  }

  # Turnover component
  if ("beta_sim" %in% names(df)) {
    f1b_t <- "beta_sim ~ rao_q_z * chm_cv_z + (1|domain/siteID)"
    c1b_t <- fit_and_report(f1b_t, df, "M1b Turnover")
    if (!is.null(c1b_t)) all_coefs[[length(all_coefs) + 1]] <- c1b_t
  }
} else {
  cat("\nSkipping Model 1b: insufficient data\n")
}


# ════════════════════════════════════════════════════════════════════════════
# Model 2: Productivity ~ RS diversity * heterogeneity
# ════════════════════════════════════════════════════════════════════════════

has_productivity <- "cumulative_mean_z" %in% names(df) && sum(!is.na(df$cumulative_mean_z)) > 10
has_functional <- "func_FRic_z" %in% names(df) && sum(!is.na(df$func_FRic_z)) > 10

if (has_productivity && has_structural) {
  f2 <- "cumulative_mean ~ chm_cv_z + (1|domain)"
  c2 <- fit_and_report(f2, df %>% distinct(siteID, .keep_all = TRUE), "M2 Productivity Basic")
  if (!is.null(c2)) all_coefs[[length(all_coefs) + 1]] <- c2

  if (has_heterogeneity) {
    f2b <- "cumulative_mean ~ chm_cv_z * Sa_500m_z + (1|domain)"
    c2b <- fit_and_report(f2b, df %>% distinct(siteID, .keep_all = TRUE),
                          "M2 Productivity x Heterogeneity")
    if (!is.null(c2b)) all_coefs[[length(all_coefs) + 1]] <- c2b
  }

  if (has_functional) {
    f2c <- "cumulative_mean ~ func_FRic_z * chm_cv_z + (1|domain)"
    c2c <- fit_and_report(f2c, df %>% distinct(siteID, .keep_all = TRUE),
                          "M2 Productivity ~ Functional Div")
    if (!is.null(c2c)) all_coefs[[length(all_coefs) + 1]] <- c2c
  }
} else {
  cat("\nSkipping Model 2: insufficient data\n")
  cat(sprintf("  has_productivity=%s, has_structural=%s\n", has_productivity, has_structural))
}


# ════════════════════════════════════════════════════════════════════════════
# Save all coefficients
# ════════════════════════════════════════════════════════════════════════════

if (length(all_coefs) > 0) {
  coefs_combined <- bind_rows(all_coefs)
  coef_path <- file.path(base_dir, "model_coefficients.csv")
  write.csv(coefs_combined, coef_path, row.names = FALSE)
  cat(sprintf("\n\nAll coefficients saved: %s (%d rows)\n", coef_path, nrow(coefs_combined)))
} else {
  cat("\nNo models were successfully fitted.\n")
}

cat("\nDone.\n")
