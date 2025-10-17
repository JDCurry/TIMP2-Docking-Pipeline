#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# ---------------- EDIT THESE ----------------
INPUT_SMI = r"C:\data\timp2\results\timp2_filtered_filtered.smi"
OUT_DIR   = r"C:\data\timp2\cns_gate"
MW_MIN, MW_MAX = 200, 550
LOGP_MIN, LOGP_MAX = 1.5, 4.0
TPSA_MAX = 90
HBD_MAX  = 2
# --------------------------------------------

Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
out_pass = Path(OUT_DIR, "timp2_CNS_pass.smi").open("w", encoding="utf-8")
out_fail = Path(OUT_DIR, "timp2_CNS_fail.smi").open("w", encoding="utf-8")

def cns_ok(m):
    mw   = Descriptors.MolWt(m)
    if not (MW_MIN <= mw <= MW_MAX): return False
    logp = Descriptors.MolLogP(m)
    tpsa = rdMolDescriptors.CalcTPSA(m)
    hbd  = rdMolDescriptors.CalcNumHBD(m)
    return (TPSA_MAX >= tpsa) and (HBD_MAX >= hbd) and (LOGP_MIN <= logp <= LOGP_MAX)

count, kept = 0, 0
with open(INPUT_SMI, "r", encoding="utf-8", errors="ignore") as fin:
    for line in fin:
        line = line.strip()
        if not line: continue
        parts = line.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else ""
        m = Chem.MolFromSmiles(smi)
        if m is None:
            out_fail.write(line + "\n")
            continue
        if cns_ok(m):
            out_pass.write(line + "\n"); kept += 1
        else:
            out_fail.write(line + "\n")
        count += 1
        if count % 1_000_000 == 0:
            print(f"Processed {count:,} lines; kept {kept:,}")

out_pass.close(); out_fail.close()
print(f"[DONE] Processed {count:,}; CNS pass {kept:,}")
print("Pass file:", Path(OUT_DIR, "timp2_CNS_pass.smi"))
