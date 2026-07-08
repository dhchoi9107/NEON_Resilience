#!/bin/bash
cd "C:/Users/star1/Documents/GitHub/NEON_Resilience"
PY="/c/Users/star1/anaconda3/python.exe"
LOG="NEON_v2/scripts/obj23/downstream.log"
echo "waiting for sentinel extraction..." > $LOG
# 추출 완료 대기 (ALL DONE 또는 19개 csv)
while true; do
  n=$(ls NEON_v2/data/sentinel/s2_*.csv 2>/dev/null | wc -l)
  if grep -q "ALL DONE" NEON_v2/scripts/obj23/extract_all2.log 2>/dev/null; then break; fi
  sleep 30
done
echo "=== extraction done ($(ls NEON_v2/data/sentinel/s2_*.csv|wc -l) sites). running downstream ===" >> $LOG
PYTHONIOENCODING=utf-8 "$PY" NEON_v2/scripts/obj23/20_compute_dhi.py >> $LOG 2>&1
echo "--- DHI done ---" >> $LOG
PYTHONIOENCODING=utf-8 "$PY" NEON_v2/scripts/obj23/40_disturbance_selfdetect.py >> $LOG 2>&1
echo "--- selfdetect done ---" >> $LOG
PYTHONIOENCODING=utf-8 "$PY" NEON_v2/scripts/obj23/50_models_obj23.py >> $LOG 2>&1
echo "--- models done ---" >> $LOG
PYTHONIOENCODING=utf-8 "$PY" NEON_v2/scripts/obj23/51_figures_obj23.py >> $LOG 2>&1
echo "=== ALL DOWNSTREAM DONE ===" >> $LOG
