#!/usr/bin/env Rscript
# ============================================================================
# Comprehensive Models: Alpha, Beta, Heterogeneity interactions
# ============================================================================

library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)

base_dir <- "E:/neon_lidar/model_results"
df <- read.csv(file.path(base_dir, "plot_level_all.csv"), stringsAsFactors = FALSE)
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)

# Z-score
zvars <- c("vert_sd_mean","FHD_mean","deepgap_fraction_mean","GC_mean",
           "rao_q","spectral_cv","spectral_FRic",
           "func_FRic","func_RaoQ",
           "chm_cv","Sa_500m","Sq_500m","cumulative_mean","top_rugosity_mean",
           "trend_mean","variation_mean")
for (v in zvars) {
  zv <- paste0(v, "_z")
  if (v %in% names(df) && sum(!is.na(df[[v]])) > 10)
    df[[zv]] <- scale(df[[v]])[, 1]
}

# Site-level aggregation
site_df <- df %>%
  group_by(siteID, domain) %>%
  summarise(across(where(is.numeric), ~mean(.x, na.rm=TRUE)), .groups="drop") %>%
  as.data.frame()
site_df$siteID <- as.factor(site_df$siteID)
site_df$domain <- as.factor(site_df$domain)
for (v in zvars) {
  zv <- paste0(v, "_z")
  if (v %in% names(site_df) && sum(!is.na(site_df[[v]])) > 5)
    site_df[[zv]] <- scale(site_df[[v]])[, 1]
}

all_coefs <- list()

fit <- function(formula_str, data, name, category) {
  cat(sprintf("  %-50s", name))
  tryCatch({
    m <- lmer(as.formula(formula_str), data = data, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
    r2 <- r.squaredGLMM(m)
    fe <- as.data.frame(coef(summary(m)))
    fe$term <- rownames(fe)
    fe$model <- name
    fe$category <- category
    fe$R2m <- r2[1,"R2m"]
    fe$R2c <- r2[1,"R2c"]
    fe$n <- nobs(m)
    cat(sprintf("  R2m=%.3f  R2c=%.3f  n=%d\n", r2[1,"R2m"], r2[1,"R2c"], nobs(m)))
    all_coefs[[length(all_coefs)+1]] <<- fe
  }, error = function(e) cat(sprintf("  ERROR: %s\n", e$message)))
}

# ═══════════════════════════════════════
cat("=== ALPHA DIVERSITY ===\n")
# ═══════════════════════════════════════

# Single predictors
fit("shannon ~ deepgap_fraction_mean_z + (1|domain/siteID)", df,
    "Deep Gap", "Alpha-Single")
fit("shannon ~ GC_mean_z + (1|domain/siteID)", df,
    "Gini", "Alpha-Single")
fit("shannon ~ FHD_mean_z + (1|domain/siteID)", df,
    "FHD", "Alpha-Single")
fit("shannon ~ vert_sd_mean_z + (1|domain/siteID)", df,
    "Vert SD", "Alpha-Single")
fit("shannon ~ rao_q_z + (1|domain/siteID)", df,
    "Spectral Rao Q", "Alpha-Single")

# Multivariate
fit("shannon ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + (1|domain/siteID)", df,
    "RaoQ + Gap + Gini", "Alpha-Multi")
fit("shannon ~ rao_q_z * FHD_mean_z + (1|domain/siteID)", df,
    "RaoQ x FHD", "Alpha-Interaction")

# With heterogeneity
fit("shannon ~ rao_q_z + GC_mean_z + chm_cv_z + (1|domain/siteID)", df,
    "RaoQ + Gini + CHM_CV", "Alpha-Heterogeneity")
fit("shannon ~ rao_q_z + GC_mean_z + Sa_500m_z + (1|domain/siteID)", df,
    "RaoQ + Gini + Sa", "Alpha-Heterogeneity")
fit("shannon ~ rao_q_z * chm_cv_z + GC_mean_z + (1|domain/siteID)", df,
    "RaoQ x CHM_CV + Gini", "Alpha-Heterogeneity")
fit("shannon ~ deepgap_fraction_mean_z * chm_cv_z + (1|domain/siteID)", df,
    "Gap x CHM_CV", "Alpha-Heterogeneity")
fit("shannon ~ deepgap_fraction_mean_z * Sa_500m_z + (1|domain/siteID)", df,
    "Gap x Sa", "Alpha-Heterogeneity")

# Full alpha model
fit("shannon ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + chm_cv_z + rao_q_z:FHD_mean_z + (1|domain/siteID)", df,
    "Full Alpha", "Alpha-Full")

# Richness
fit("richness ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + chm_cv_z + (1|domain/siteID)", df,
    "Richness Full", "Richness")


# ═══════════════════════════════════════
cat("\n=== BETA DIVERSITY (site-level) ===\n")
# ═══════════════════════════════════════

# Bray-Curtis
fit("bray_mean ~ rao_q_z + (1|domain)", site_df,
    "Bray ~ RaoQ", "Beta-Bray")
fit("bray_mean ~ deepgap_fraction_mean_z + (1|domain)", site_df,
    "Bray ~ Gap", "Beta-Bray")
fit("bray_mean ~ GC_mean_z + (1|domain)", site_df,
    "Bray ~ Gini", "Beta-Bray")
fit("bray_mean ~ chm_cv_z + (1|domain)", site_df,
    "Bray ~ CHM_CV", "Beta-Bray")
fit("bray_mean ~ rao_q_z + GC_mean_z + chm_cv_z + (1|domain)", site_df,
    "Bray ~ RaoQ + Gini + CHM_CV", "Beta-Bray-Multi")
fit("bray_mean ~ rao_q_z * chm_cv_z + (1|domain)", site_df,
    "Bray ~ RaoQ x CHM_CV", "Beta-Bray-Het")

# Turnover
fit("beta_sim ~ rao_q_z + (1|domain)", site_df,
    "Turnover ~ RaoQ", "Beta-Turnover")
fit("beta_sim ~ deepgap_fraction_mean_z + (1|domain)", site_df,
    "Turnover ~ Gap", "Beta-Turnover")
fit("beta_sim ~ GC_mean_z + (1|domain)", site_df,
    "Turnover ~ Gini", "Beta-Turnover")
fit("beta_sim ~ rao_q_z + GC_mean_z + (1|domain)", site_df,
    "Turnover ~ RaoQ + Gini", "Beta-Turnover-Multi")
fit("beta_sim ~ rao_q_z * chm_cv_z + (1|domain)", site_df,
    "Turnover ~ RaoQ x CHM_CV", "Beta-Turnover-Het")

# Nestedness
fit("beta_sne ~ rao_q_z + (1|domain)", site_df,
    "Nestedness ~ RaoQ", "Beta-Nestedness")
fit("beta_sne ~ deepgap_fraction_mean_z + (1|domain)", site_df,
    "Nestedness ~ Gap", "Beta-Nestedness")
fit("beta_sne ~ GC_mean_z + (1|domain)", site_df,
    "Nestedness ~ Gini", "Beta-Nestedness")


# ═══════════════════════════════════════
cat("\n=== PRODUCTIVITY (site-level) ===\n")
# ═══════════════════════════════════════

fit("cumulative_mean ~ top_rugosity_mean_z + (1|domain)", site_df,
    "Cumul ~ Rugosity", "Productivity")
fit("cumulative_mean ~ FHD_mean_z + (1|domain)", site_df,
    "Cumul ~ FHD", "Productivity")
fit("cumulative_mean ~ rao_q_z + (1|domain)", site_df,
    "Cumul ~ RaoQ", "Productivity")
fit("cumulative_mean ~ top_rugosity_mean_z + chm_cv_z + (1|domain)", site_df,
    "Cumul ~ Rugosity + CHM_CV", "Productivity-Het")
fit("cumulative_mean ~ top_rugosity_mean_z * Sa_500m_z + (1|domain)", site_df,
    "Cumul ~ Rugosity x Sa", "Productivity-Het")
fit("trend_mean ~ top_rugosity_mean_z + (1|domain)", site_df,
    "Trend ~ Rugosity", "Productivity-Trend")
fit("trend_mean ~ chm_cv_z + (1|domain)", site_df,
    "Trend ~ CHM_CV", "Productivity-Trend")
fit("variation_mean ~ chm_cv_z + (1|domain)", site_df,
    "Variation ~ CHM_CV", "Productivity-Stability")


# ═══════════════════════════════════════
# Save
# ═══════════════════════════════════════
if (length(all_coefs) > 0) {
  coefs <- bind_rows(all_coefs)
  write.csv(coefs, file.path(base_dir, "comprehensive_model_coefficients.csv"), row.names = FALSE)
  cat(sprintf("\nSaved: %d coefficients from %d models\n",
              nrow(coefs), length(unique(coefs$model))))
}
cat("Done.\n")
