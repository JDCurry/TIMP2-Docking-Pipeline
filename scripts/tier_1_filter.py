# -*- coding: utf-8 -*-
"""
Strict Tier 1 filtering for ADMET-AI outputs.
Usage:
    python strict_tier_filter.py input_file.csv
Outputs:
    input_file_tiered.csv
    input_file_tier1_only.csv
"""

import sys
import pandas as pd
import numpy as np

# If no argument provided, fallback to a test filename
if len(sys.argv) < 2:
    infile = "preds.csv"
else:
    infile = sys.argv[1]

df = pd.read_csv(infile)

# --- Helper: boolean casting for numeric or binary columns ---
def flag(x):
    try:
        return int(float(x) > 0.5)
    except:
        return 0

# --- Define hard and soft flags ---
hard_flags = ["AMES", "hERG", "Carcinogens_Lagunin", "DILI"]
soft_flags = ["ClinTox", "SR-ARE", "SR-HSE", "SR-p53", "NR-AhR"]

# --- Build composite penalty ---
df["hard_penalty"] = df[hard_flags].applymap(flag).sum(axis=1)
df["soft_penalty"] = df[soft_flags].applymap(flag).sum(axis=1)

# --- Assign safety tiers ---
df["tier"] = np.select(
    [
        (df["hard_penalty"] == 0)
        & (df["soft_penalty"] == 0)
        & (df["ClinTox"] < 0.6)
        & (df["LD50_Zhu"] > 2000),
        (df["hard_penalty"] == 0)
        & ((df["soft_penalty"] <= 2) | (df["ClinTox"] < 0.7)),
        (df["hard_penalty"] > 0) | (df["ClinTox"] >= 0.7) | (df["LD50_Zhu"] <= 2000),
    ],
    [1, 2, 3],
)

# --- Export results ---
out_tiered = infile.replace(".csv", "_tiered.csv")
out_tier1 = infile.replace(".csv", "_tier1_only.csv")

df.to_csv(out_tiered, index=False, encoding="utf-8-sig")
df[df["tier"] == 1].to_csv(out_tier1, index=False, encoding="utf-8-sig")

# --- Print summary ---
summary = df["tier"].value_counts().sort_index()
print("Tier summary:")
for t, n in summary.items():
    label = {1: "Tier 1 – clean", 2: "Tier 2 – caution", 3: "Tier 3 – fail"}[t]
    print(f"  {t}: {n} compounds ({label})")

print(f"\nWrote:\n  {out_tiered}\n  {out_tier1}")
