#!/bin/bash
export MSYS_NO_PATHCONV=1
SRC='E:\neon_lidar'
DST='\10.10.170.55\home\dhchoi_main\NEON_Resilience_data'
LOG="scripts_pipeline/_pipeline_state/copy_neon.log"
echo "copy start $(date)" > "$LOG"
for F in structural_diversity structural_diversity_1m structural_diversity_1m_plots stand_age spectral_diversity hyperspectral_brdf_corrected hyperspectral_plots_001 rgb_plots vegetation_structure env_heterogeneity functional_diversity taxonomic_diversity model_results productivity_dhi; do
  echo "===== $F $(date +%H:%M:%S) =====" >> "$LOG"
  robocopy "$SRC\$F" "$DST\$F" /E /MT:16 /R:1 /W:3 /NFL /NDL /NJH /NP /BYTES >> "$LOG" 2>&1
  echo "  rc=$?" >> "$LOG"
done
echo "copy done $(date)" >> "$LOG"
