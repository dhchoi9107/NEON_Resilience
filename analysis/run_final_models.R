#!/usr/bin/env Rscript
# ============================================================================
# Final Models with VIF-checked variable set
# ============================================================================
# Variables retained (all VIF < 5):
#   Structural (4): vert_sd, FHD, deepgap_fraction, GC
#   Spectral   (3): rao_q, spectral_cv, spectral_FRic
#   Functional (2): func_FRic, func_RaoQ
# ============================================================================

library(lme4)
library(lmerTest)
library(MuMIn)
library(dplyr)

base_dir <- "E:/neon_lidar/model_results"
df <- read.csv(file.path(base_dir, "plot_level_all.csv"), stringsAsFactors = FALSE)
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)
cat(sprintf("Data: %d plots, %d sites\n\n", nrow(df), length(unique(df$siteID))))

# Z-score
zvars <- c("vert_sd_mean","FHD_mean","deepgap_fraction_mean","GC_mean",
           "rao_q","spectral_cv","spectral_FRic",
           "func_FRic","func_RaoQ",
           "chm_cv","Sa_500m","cumulative_mean","top_rugosity_mean")
for (v in zvars) {
  zv <- paste0(v, "_z")
  if (v %in% names(df) && sum(!is.na(df[[v]])) > 10) {
    df[[zv]] <- scale(df[[v]])[, 1]
  }
}

all_coefs <- list()

fit <- function(formula_str, data, name) {
  cat(sprintf("── %s ──\n", name))
  tryCatch({
    m <- lmer(as.formula(formula_str), data = data, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
    r2 <- r.squaredGLMM(m)
    fe <- as.data.frame(coef(summary(m)))
    fe$term <- rownames(fe)
    fe$model <- name
    fe$R2m <- r2[1,"R2m"]
    fe$R2c <- r2[1,"R2c"]
    fe$n <- nobs(m)

    # Print
    cat(sprintf("  R2m=%.3f  R2c=%.3f  n=%d\n", r2[1,"R2m"], r2[1,"R2c"], nobs(m)))
    for (i in 1:nrow(fe)) {
      p <- fe[["Pr(>|t|)"]][i]
      s <- ifelse(p<0.001,"***", ifelse(p<0.01,"**", ifelse(p<0.05,"*", ifelse(p<0.1,".",""))))
      cat(sprintf("    %-35s  b=%+7.4f  p=%.4f %s\n", fe$term[i], fe$Estimate[i], p, s))
    }
    cat("\n")
    all_coefs[[length(all_coefs)+1]] <<- fe
  }, error = function(e) cat(sprintf("  ERROR: %s\n\n", e$message)))
}

# ═══════════════════════════════════════
cat("=== ALPHA DIVERSITY MODELS ===\n\n")
# ═══════════════════════════════════════

# Individual structural predictors
fit("shannon ~ vert_sd_mean_z + (1|domain/siteID)", df,
    "Alpha ~ Vert SD")
fit("shannon ~ FHD_mean_z + (1|domain/siteID)", df,
    "Alpha ~ FHD")
fit("shannon ~ deepgap_fraction_mean_z + (1|domain/siteID)", df,
    "Alpha ~ Deep Gap")
fit("shannon ~ GC_mean_z + (1|domain/siteID)", df,
    "Alpha ~ Gini")

# Multivariate structural (no collinear pairs)
fit("shannon ~ FHD_mean_z + deepgap_fraction_mean_z + GC_mean_z + (1|domain/siteID)", df,
    "Alpha ~ FHD + Gap + Gini")

# Spectral
fit("shannon ~ rao_q_z + (1|domain/siteID)", df,
    "Alpha ~ Rao Q")

# Spectral x Structural interaction
fit("shannon ~ rao_q_z * FHD_mean_z + (1|domain/siteID)", df,
    "Alpha ~ RaoQ * FHD")
fit("shannon ~ rao_q_z * deepgap_fraction_mean_z + (1|domain/siteID)", df,
    "Alpha ~ RaoQ * Gap")

# Full model: spectral + structural + interaction
fit("shannon ~ rao_q_z + FHD_mean_z + deepgap_fraction_mean_z + GC_mean_z + rao_q_z:FHD_mean_z + (1|domain/siteID)", df,
    "Alpha ~ Full (Spectral + Structural)")

# Functional
fit("shannon ~ func_FRic_z + func_RaoQ_z + (1|domain/siteID)", df,
    "Alpha ~ Func FRic + RaoQ")


# ═══════════════════════════════════════
cat("\n=== RICHNESS MODELS ===\n\n")
# ═══════════════════════════════════════

fit("richness ~ rao_q_z * FHD_mean_z + deepgap_fraction_mean_z + GC_mean_z + (1|domain/siteID)", df,
    "Richness ~ Full")


# ═══════════════════════════════════════
cat("\n=== BETA DIVERSITY MODELS (site-level) ===\n\n")
# ═══════════════════════════════════════

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

fit("bray_mean ~ rao_q_z + FHD_mean_z + deepgap_fraction_mean_z + (1|domain)", site_df,
    "Beta Bray ~ Spectral + Structural")
fit("beta_sim ~ rao_q_z + FHD_mean_z + deepgap_fraction_mean_z + (1|domain)", site_df,
    "Beta Turnover ~ Spectral + Structural")


# ═══════════════════════════════════════
cat("\n=== PRODUCTIVITY MODELS (site-level) ===\n\n")
# ═══════════════════════════════════════

fit("cumulative_mean ~ FHD_mean_z + deepgap_fraction_mean_z + (1|domain)", site_df,
    "Productivity ~ FHD + Gap")
fit("cumulative_mean ~ rao_q_z + FHD_mean_z + Sa_500m_z + (1|domain)", site_df,
    "Productivity ~ RaoQ + FHD + Sa")
fit("cumulative_mean ~ top_rugosity_mean_z + (1|domain)", site_df,
    "Productivity ~ Rugosity")


# ═══════════════════════════════════════
# Save
# ═══════════════════════════════════════
if (length(all_coefs) > 0) {
  coefs <- bind_rows(all_coefs)
  write.csv(coefs, file.path(base_dir, "final_model_coefficients.csv"), row.names = FALSE)
  cat(sprintf("\nSaved: %d coefficients from %d models\n",
              nrow(coefs), length(unique(coefs$model))))
}
cat("Done.\n")
