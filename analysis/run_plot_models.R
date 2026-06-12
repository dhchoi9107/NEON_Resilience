#!/usr/bin/env Rscript
# ============================================================================
# Plot-level Mixed-Effects Models (n=637 plots)
# ============================================================================
# Tests all diversity relationships at plot level with 1m resolution data.
# This gives much more statistical power than site-level (n~19) models.
#
# Models:
#   1) Alpha ~ structural diversity + (1|domain/site)
#   2) Alpha ~ spectral diversity + (1|domain/site)
#   3) Alpha ~ spectral * structural + (1|domain/site)
#   4) Alpha ~ functional diversity + (1|domain/site)
#   5) Beta ~ structural + spectral + (1|domain)
#   6) Productivity ~ diversity * heterogeneity + (1|domain)
# ============================================================================

library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(dplyr)

base_dir <- "E:/neon_lidar/model_results"
fig_dir <- file.path(base_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load data ──
cat("Loading plot-level data...\n")
df <- read.csv(file.path(base_dir, "plot_level_all.csv"), stringsAsFactors = FALSE)
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)
cat(sprintf("  %d plots, %d sites, %d domains\n", nrow(df), length(unique(df$siteID)), length(unique(df$domain))))

# ── Z-score all predictors ──
pred_cols <- c("top_rugosity_mean", "mean_max_canopy_ht_mean", "vert_sd_mean",
               "vertCV_mean", "FHD_mean", "LAI_mean", "GC_mean", "deepgap_fraction_mean",
               "max_canopy_ht_mean", "HeightRatio_mean",
               "rao_q", "spectral_cv", "spectral_shannon",
               "spectral_FRic", "spectral_FDiv", "spectral_FEve",
               "func_FRic", "func_FDiv", "func_FEve", "func_RaoQ",
               "chm_cv", "Sa_500m", "Sq_500m",
               "cumulative_mean", "trend_mean", "variation_mean")

for (col in pred_cols) {
  zcol <- paste0(col, "_z")
  if (col %in% names(df) && sum(!is.na(df[[col]])) > 10) {
    df[[zcol]] <- scale(df[[col]])[, 1]
  }
}

# ── Fit function ──
all_coefs <- list()

fit_lmer <- function(formula_str, data, model_name) {
  cat(sprintf("\n══ %s ══\n", model_name))
  cat(sprintf("  %s\n", formula_str))

  tryCatch({
    m <- lmer(as.formula(formula_str), data = data, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))

    s <- summary(m)
    r2 <- r.squaredGLMM(m)

    cat(sprintf("  R2m=%.3f  R2c=%.3f  n=%d\n", r2[1,"R2m"], r2[1,"R2c"], nobs(m)))

    # Print fixed effects
    fe <- as.data.frame(s$coefficients)
    fe$sig <- ifelse(fe[["Pr(>|t|)"]] < 0.001, "***",
              ifelse(fe[["Pr(>|t|)"]] < 0.01, "**",
              ifelse(fe[["Pr(>|t|)"]] < 0.05, "*",
              ifelse(fe[["Pr(>|t|)"]] < 0.1, ".", ""))))
    for (i in 1:nrow(fe)) {
      cat(sprintf("    %-35s  b=%+7.3f  SE=%.3f  p=%.4f %s\n",
          rownames(fe)[i], fe$Estimate[i], fe$`Std. Error`[i], fe$`Pr(>|t|)`[i], fe$sig[i]))
    }

    # Save coefficients
    coef_df <- as.data.frame(s$coefficients)
    coef_df$term <- rownames(coef_df)
    coef_df$model <- model_name
    coef_df$R2m <- r2[1,"R2m"]
    coef_df$R2c <- r2[1,"R2c"]
    coef_df$n <- nobs(m)
    all_coefs[[length(all_coefs) + 1]] <<- coef_df

    return(m)
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
    return(NULL)
  })
}


# ════════════════════════════════════════════════════════════════
# MODEL 1: Alpha ~ Structural diversity
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 1: Alpha ~ Structural\n", strrep("=", 60), "\n")

fit_lmer("shannon ~ top_rugosity_mean_z + (1|domain/siteID)", df, "Alpha ~ Top Rugosity")
fit_lmer("shannon ~ vert_sd_mean_z + (1|domain/siteID)", df, "Alpha ~ Vertical SD")
fit_lmer("shannon ~ FHD_mean_z + (1|domain/siteID)", df, "Alpha ~ FHD")
fit_lmer("shannon ~ deepgap_fraction_mean_z + (1|domain/siteID)", df, "Alpha ~ Deep Gap Frac")
fit_lmer("shannon ~ GC_mean_z + (1|domain/siteID)", df, "Alpha ~ Gini")
fit_lmer("shannon ~ mean_max_canopy_ht_mean_z + (1|domain/siteID)", df, "Alpha ~ Mean Canopy Ht")
fit_lmer("shannon ~ top_rugosity_mean_z + FHD_mean_z + deepgap_fraction_mean_z + (1|domain/siteID)", df,
         "Alpha ~ Rugosity + FHD + Gap")


# ════════════════════════════════════════════════════════════════
# MODEL 2: Alpha ~ Spectral diversity
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 2: Alpha ~ Spectral\n", strrep("=", 60), "\n")

fit_lmer("shannon ~ rao_q_z + (1|domain/siteID)", df, "Alpha ~ Spectral Rao Q")
fit_lmer("shannon ~ spectral_FRic_z + (1|domain/siteID)", df, "Alpha ~ Spectral FRic")
fit_lmer("shannon ~ spectral_shannon_z + (1|domain/siteID)", df, "Alpha ~ Spectral Shannon")


# ════════════════════════════════════════════════════════════════
# MODEL 3: Alpha ~ Spectral * Structural (interaction)
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 3: Alpha ~ Spectral x Structural\n", strrep("=", 60), "\n")

fit_lmer("shannon ~ rao_q_z * top_rugosity_mean_z + (1|domain/siteID)", df,
         "Alpha ~ RaoQ x Rugosity")
fit_lmer("shannon ~ rao_q_z * FHD_mean_z + (1|domain/siteID)", df,
         "Alpha ~ RaoQ x FHD")
fit_lmer("shannon ~ spectral_FRic_z * vert_sd_mean_z + (1|domain/siteID)", df,
         "Alpha ~ SpecFRic x VertSD")


# ════════════════════════════════════════════════════════════════
# MODEL 4: Alpha ~ Functional diversity (combined trait space)
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 4: Alpha ~ Functional Diversity\n", strrep("=", 60), "\n")

fit_lmer("shannon ~ func_FRic_z + (1|domain/siteID)", df, "Alpha ~ Func FRic")
fit_lmer("shannon ~ func_RaoQ_z + (1|domain/siteID)", df, "Alpha ~ Func Rao Q")
fit_lmer("shannon ~ func_FRic_z + func_FDiv_z + (1|domain/siteID)", df, "Alpha ~ Func FRic + FDiv")


# ════════════════════════════════════════════════════════════════
# MODEL 5: Beta diversity (site-level in plot data)
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 5: Beta ~ RS diversity\n", strrep("=", 60), "\n")

# Use site-level means
site_df <- df %>%
  group_by(siteID, domain) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm=TRUE)), .groups="drop") %>%
  as.data.frame()
site_df$siteID <- as.factor(site_df$siteID)
site_df$domain <- as.factor(site_df$domain)

# Re-scale at site level
for (col in pred_cols) {
  zcol <- paste0(col, "_z")
  if (col %in% names(site_df) && sum(!is.na(site_df[[col]])) > 5) {
    site_df[[zcol]] <- scale(site_df[[col]])[, 1]
  }
}

fit_lmer("bray_mean ~ rao_q_z + top_rugosity_mean_z + (1|domain)", site_df,
         "Beta Bray ~ RaoQ + Rugosity (site)")
fit_lmer("beta_sim ~ rao_q_z + top_rugosity_mean_z + (1|domain)", site_df,
         "Beta Turnover ~ RaoQ + Rugosity (site)")
fit_lmer("bray_mean ~ func_FRic_z + chm_cv_z + (1|domain)", site_df,
         "Beta Bray ~ FuncFRic + CHM_CV (site)")


# ════════════════════════════════════════════════════════════════
# MODEL 6: Productivity ~ diversity * heterogeneity
# ════════════════════════════════════════════════════════════════
cat("\n", strrep("=", 60), "\n  MODEL GROUP 6: Productivity\n", strrep("=", 60), "\n")

fit_lmer("cumulative_mean ~ top_rugosity_mean_z + (1|domain)", site_df,
         "Productivity ~ Rugosity (site)")
fit_lmer("cumulative_mean ~ top_rugosity_mean_z * Sa_500m_z + (1|domain)", site_df,
         "Productivity ~ Rugosity x Sa (site)")
fit_lmer("cumulative_mean ~ func_FRic_z + Sa_500m_z + (1|domain)", site_df,
         "Productivity ~ FuncFRic + Sa (site)")


# ════════════════════════════════════════════════════════════════
# Save all coefficients
# ════════════════════════════════════════════════════════════════
if (length(all_coefs) > 0) {
  coefs <- bind_rows(all_coefs)
  coef_path <- file.path(base_dir, "plot_model_coefficients.csv")
  write.csv(coefs, coef_path, row.names = FALSE)
  cat(sprintf("\n\nAll coefficients: %s (%d rows)\n", coef_path, nrow(coefs)))
}

cat("\nDone.\n")
