# split_smi.py
from pathlib import Path
INPUT  = r"C:\data\timp2\cns_gate\CNS_STRICT.smi"
OUTDIR = r"C:\data\timp2\chunks_strict"
CHUNK  = 100_000
Path(OUTDIR).mkdir(parents=True, exist_ok=True)
i = 0; n = 0; w = None
with open(INPUT, "r", encoding="utf-8", errors="ignore") as fin:
    for line in fin:
        if n % CHUNK == 0:
            if w: w.close()
            i += 1
            w = open(Path(OUTDIR)/f"strict_{i:04d}.smi", "w", encoding="utf-8")
        w.write(line); n += 1
if w: w.close()
print(f"Chunks: {i}, size≈{CHUNK}")
