library(lme4); library(lmerTest); library(MuMIn); library(dplyr)

df <- read.csv("E:/neon_lidar/model_results/plot_level_complete.csv")
df$siteID <- as.factor(df$siteID)
df$domain <- as.factor(df$domain)

zvars <- c("deepgap_fraction_mean","GC_mean","FHD_mean","rao_q","vert_sd_mean",
           "deepgap_fraction_sd","mean_max_canopy_ht_sd")
for (v in zvars) df[[paste0(v,"_z")]] <- scale(df[[v]])[,1]

all_coefs <- list()
fit <- function(f, data, name, cat_label) {
  tryCatch({
    m <- lmer(as.formula(f), data=data, REML=TRUE, control=lmerControl(optimizer="bobyqa"))
    r2 <- r.squaredGLMM(m)
    fe <- as.data.frame(coef(summary(m)))
    fe$term <- rownames(fe); fe$model <- name; fe$category <- cat_label
    fe$R2m <- r2[1,"R2m"]; fe$R2c <- r2[1,"R2c"]; fe$n <- nobs(m)
    fe$AIC <- AIC(m)
    cat(sprintf("%-50s R2m=%.3f R2c=%.3f AIC=%.0f\n", name, r2[1,"R2m"], r2[1,"R2c"], AIC(m)))
    for (i in 1:nrow(fe)) {
      if (fe$term[i]=="(Intercept)") next
      p <- fe[["Pr(>|t|)"]][i]
      s <- ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".",""))))
      cat(sprintf("  %-40s b=%+7.4f p=%.4f %s\n", fe$term[i], fe$Estimate[i], p, s))
    }
    all_coefs[[length(all_coefs)+1]] <<- fe
  }, error=function(e) cat(sprintf("%-50s ERROR: %s\n", name, e$message)))
}

cat("=== RANDOM INTERCEPT vs RANDOM SLOPE ===\n\n")

fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + rao_q_z + (1|domain/siteID)", df,
    "M1: Intercept only", "compare")

fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + rao_q_z + (1+deepgap_fraction_mean_z|siteID)", df,
    "M2: Slope on gap", "compare")

fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + rao_q_z + (1+deepgap_fraction_mean_z+GC_mean_z||siteID)", df,
    "M3: Slopes on gap+GC", "compare")

# Full model with heterogeneity + random slopes
fit("shannon ~ deepgap_fraction_mean_z + GC_mean_z + rao_q_z + deepgap_fraction_sd_z + (1+deepgap_fraction_mean_z||siteID)", df,
    "M4: M2 + Gap_sd", "compare")

cat("\n=== WITHIN-SITE ANALYSIS (site-centered) ===\n\n")

# Center within site to isolate within-site relationships
df2 <- df %>% group_by(siteID) %>% mutate(
  gap_wc = as.numeric(scale(deepgap_fraction_mean)),
  gc_wc  = as.numeric(scale(GC_mean)),
  raoq_wc = as.numeric(scale(rao_q)),
  fhd_wc = as.numeric(scale(FHD_mean)),
  gapsd_wc = as.numeric(scale(deepgap_fraction_sd)),
  shannon_wc = as.numeric(scale(shannon))
) %>% ungroup() %>% as.data.frame()
df2 <- df2[complete.cases(df2[,c("shannon","gap_wc","gc_wc","raoq_wc")]),]

m_wc <- lm(shannon ~ gap_wc + gc_wc + raoq_wc, df2)
cat(sprintf("Within-site OLS: R2=%.3f adj.R2=%.3f n=%d\n",
    summary(m_wc)$r.squared, summary(m_wc)$adj.r.squared, nrow(df2)))
fe_wc <- summary(m_wc)$coefficients
for (i in 1:nrow(fe_wc)) {
  p <- fe_wc[i,4]
  s <- ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".",""))))
  cat(sprintf("  %-30s b=%+7.4f p=%.4f %s\n", rownames(fe_wc)[i], fe_wc[i,1], p, s))
}

m_wc2 <- lm(shannon ~ gap_wc + gc_wc + raoq_wc + fhd_wc + gapsd_wc, df2)
cat(sprintf("\nWithin-site + het: R2=%.3f adj.R2=%.3f\n",
    summary(m_wc2)$r.squared, summary(m_wc2)$adj.r.squared))
fe_wc2 <- summary(m_wc2)$coefficients
for (i in 1:nrow(fe_wc2)) {
  p <- fe_wc2[i,4]
  s <- ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*",ifelse(p<0.1,".",""))))
  cat(sprintf("  %-30s b=%+7.4f p=%.4f %s\n", rownames(fe_wc2)[i], fe_wc2[i,1], p, s))
}

cat("\n=== BETA (LCBD) WITH RANDOM SLOPES ===\n\n")

fit("lcbd_bray ~ deepgap_fraction_mean_z + GC_mean_z + rao_q_z + deepgap_fraction_sd_z + (1+deepgap_fraction_mean_z||siteID)", df,
    "LCBD: slopes + het", "beta")

# Save
coefs <- bind_rows(all_coefs)
write.csv(coefs, "E:/neon_lidar/model_results/improved_model_coefficients.csv", row.names=FALSE)
cat(sprintf("\nSaved: %d coefficients\n", nrow(coefs)))
