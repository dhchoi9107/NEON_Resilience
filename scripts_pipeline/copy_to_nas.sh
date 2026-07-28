#!/bin/bash
DST="//10.10.170.55/home/dhchoi_main/NEON_Resilience_data"
LOG="scripts_pipeline/_pipeline_state/copy_nas.log"
echo "copy E: -> NAS start $(date)" > "$LOG"
for F in structural_diversity structural_diversity_1m structural_diversity_1m_plots stand_age spectral_diversity hyperspectral_brdf_corrected hyperspectral_plots_001 rgb_plots env_heterogeneity functional_diversity taxonomic_diversity model_results productivity_dhi; do
  echo "===== $F  start $(date +%H:%M:%S) =====" >> "$LOG"
  cp -rn "/e/neon_lidar/$F" "$DST/" 2>>"$LOG"
  echo "  done rc=$? $(date +%H:%M:%S)" >> "$LOG"
done
echo "copy done $(date)" >> "$LOG"
