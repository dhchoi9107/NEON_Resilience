# -*- coding: utf-8 -*-
from faster_whisper import WhisperModel
src=r"C:\Users\star1\Downloads\안서동 2.m4a"
out=r"C:\Users\star1\Downloads\안서동 2_transcript.txt"
# large-v3: best for code-switched Korean/English technical speech
m=WhisperModel("large-v3", device="cpu", compute_type="int8")
segs, info = m.transcribe(src, language="ko", task="transcribe", beam_size=5,
                          vad_filter=True, condition_on_previous_text=True,
                          initial_prompt="원격탐사 논문 리뷰 논의. Remote Sensing of Environment, LiDAR, GPP, beta diversity, nestedness, turnover, Bayesian, citation, supplementary 등 영어 용어가 섞여 있음.")
lines=[]
for s in segs:
    lines.append(f"[{int(s.start//60):02d}:{s.start%60:05.2f}] {s.text.strip()}")
open(out,"w",encoding="utf-8").write("\n".join(lines))
print(f"DURATION {info.duration:.1f}s SEGMENTS {len(lines)} SAVED")
