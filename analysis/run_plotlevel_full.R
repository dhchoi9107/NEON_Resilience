#!/usr/bin/env Rscript
# ============================================================================
# Plot-level models: Alpha, Beta (LCBD), Heterogeneity - all at plot level
# ============================================================================
# Beta = LCBD (Local Contribution to Beta Diversity) per plot
# Heterogeneity = within-plot spatial SD of structural metrics at 1m
# ============================================================================

library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)

base_dir <- "E:/neon_lidar/model_results"
df <- read.csv(file.path(base_dir, "plot_level_complete.csv"), stringsAsFactors = FALSE)
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)
cat(sprintf("Data: %d plots, %d sites\n\n", nrow(df), length(unique(df$siteID))))

# Z-score all predictors
zvars <- c("vert_sd_mean","FHD_mean","deepgap_fraction_mean","GC_mean",
           "rao_q","spectral_cv","spectral_FRic","spectral_FDiv",
           "func_FRic","func_RaoQ",
           # Plot-level heterogeneity (within-plot spatial SD at 1m)
           "mean_max_canopy_ht_sd","FHD_sd","vert_sd_sd","deepgap_fraction_sd",
           "top_rugosity_sd","GC_sd","LAI_sd",
           # LCBD
           "lcbd_bray")
for (v in zvars) {
  zv <- paste0(v, "_z")
  if (v %in% names(df) && sum(!is.na(df[[v]])) > 10)
    df[[zv]] <- scale(df[[v]])[, 1]
}

all_coefs <- list()
fit <- function(formula_str, data, name, category) {
  cat(sprintf("  %-55s", name))
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
    cat(sprintf("R2m=%.3f  R2c=%.3f  n=%d\n", r2[1,"R2m"], r2[1,"R2c"], nobs(m)))
    for (i in 1:nrow(fe)) {
      p <- fe[["Pr(>|t|)"]][i]
      s <- ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*", ifelse(p<0.1,".",""))))
      if (fe$term[i] != "(Intercept)")
        cat(sprintf("      %-40s b=%+7.4f  p=%.4f %s\n", fe$term[i], fe$Estimate[i], p, s))
    }
    all_coefs[[length(all_coefs)+1]] <<- fe
  }, error = function(e) cat(sprintf("ERROR: %s\n", e$message)))
}


# ═══════════════════════════════════════
cat("===== ALPHA ~ Structural + Heterogeneity =====\n")
# ═══════════════════════════════════════

# Structural diversity (mean) = composition of structure
# Structural heterogeneity (sd) = within-plot spatial variation of structure

fit("shannon ~ deepgap_fraction_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Alpha ~ Gap_mean + Gap_sd", "Alpha-Het")
fit("shannon ~ GC_mean_z + GC_sd_z + (1|domain/siteID)", df,
    "Alpha ~ Gini_mean + Gini_sd", "Alpha-Het")
fit("shannon ~ FHD_mean_z + FHD_sd_z + (1|domain/siteID)", df,
    "Alpha ~ FHD_mean + FHD_sd", "Alpha-Het")
fit("shannon ~ mean_max_canopy_ht_sd_z + (1|domain/siteID)", df,
    "Alpha ~ CHM_sd (within-plot het)", "Alpha-Het")
fit("shannon ~ vert_sd_sd_z + (1|domain/siteID)", df,
    "Alpha ~ VertSD_sd (within-plot het)", "Alpha-Het")

# Combined: structural mean + within-plot heterogeneity
fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Alpha ~ Gap + Gini + Gap_sd", "Alpha-Het-Multi")
fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + mean_max_canopy_ht_sd_z + (1|domain/siteID)", df,
    "Alpha ~ Gap + Gini + CHM_sd", "Alpha-Het-Multi")

# Spectral + heterogeneity
fit("shannon ~ rao_q_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Alpha ~ RaoQ + Gap_sd", "Alpha-Spec-Het")
fit("shannon ~ rao_q_z * deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Alpha ~ RaoQ x Gap_sd", "Alpha-Spec-Het")
fit("shannon ~ rao_q_z * mean_max_canopy_ht_sd_z + (1|domain/siteID)", df,
    "Alpha ~ RaoQ x CHM_sd", "Alpha-Spec-Het")

# Full alpha: structural + spectral + heterogeneity
fit("shannon ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Alpha FULL: RaoQ + Gap + Gini + Gap_sd", "Alpha-Full")


# ═══════════════════════════════════════
cat("\n===== BETA (LCBD) ~ Structural + Spectral + Heterogeneity =====\n")
# ═══════════════════════════════════════

# LCBD = how unique each plot is (mean Bray-Curtis to all other plots at site)
fit("lcbd_bray ~ deepgap_fraction_mean_z + (1|domain/siteID)", df,
    "LCBD ~ Gap_mean", "Beta-LCBD")
fit("lcbd_bray ~ GC_mean_z + (1|domain/siteID)", df,
    "LCBD ~ Gini_mean", "Beta-LCBD")
fit("lcbd_bray ~ FHD_mean_z + (1|domain/siteID)", df,
    "LCBD ~ FHD_mean", "Beta-LCBD")
fit("lcbd_bray ~ rao_q_z + (1|domain/siteID)", df,
    "LCBD ~ RaoQ", "Beta-LCBD")

# LCBD + heterogeneity
fit("lcbd_bray ~ deepgap_fraction_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "LCBD ~ Gap_mean + Gap_sd", "Beta-LCBD-Het")
fit("lcbd_bray ~ GC_mean_z + mean_max_canopy_ht_sd_z + (1|domain/siteID)", df,
    "LCBD ~ Gini + CHM_sd", "Beta-LCBD-Het")
fit("lcbd_bray ~ rao_q_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "LCBD ~ RaoQ + Gap_sd", "Beta-LCBD-Het")

# LCBD interactions
fit("lcbd_bray ~ rao_q_z * deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "LCBD ~ RaoQ x Gap_sd", "Beta-LCBD-Int")
fit("lcbd_bray ~ deepgap_fraction_mean_z * mean_max_canopy_ht_sd_z + (1|domain/siteID)", df,
    "LCBD ~ Gap x CHM_sd", "Beta-LCBD-Int")

# Full LCBD model
fit("lcbd_bray ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "LCBD FULL: RaoQ + Gap + Gini + Gap_sd", "Beta-LCBD-Full")


# ═══════════════════════════════════════
cat("\n===== RICHNESS =====\n")
# ═══════════════════════════════════════

fit("richness ~ rao_q_z + deepgap_fraction_mean_z + GC_mean_z + deepgap_fraction_sd_z + (1|domain/siteID)", df,
    "Richness FULL: RaoQ + Gap + Gini + Gap_sd", "Richness")


# ═══════════════════════════════════════
# Save
# ═══════════════════════════════════════
if (length(all_coefs) > 0) {
  coefs <- bind_rows(all_coefs)
  write.csv(coefs, file.path(base_dir, "plotlevel_full_coefficients.csv"), row.names = FALSE)
  cat(sprintf("\nSaved: %d coefficients from %d models\n",
              nrow(coefs), length(unique(coefs$model))))
}
cat("Done.\n")
