#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

# -------- EDIT THESE --------
INPUT_SMI = r"C:\data\timp2\results\timp2_filtered_filtered.smi"
OUT_DIR   = r"C:\data\timp2\cns_gate"
# ----------------------------

Path(OUT_DIR).mkdir(parents=True, exist_ok=True)
f_strict   = Path(OUT_DIR, "CNS_STRICT.smi").open("w", encoding="utf-8")
f_moderate = Path(OUT_DIR, "CNS_MODERATE.smi").open("w", encoding="utf-8")
f_lenient  = Path(OUT_DIR, "CNS_LENIENT.smi").open("w", encoding="utf-8")
f_fail     = Path(OUT_DIR, "CNS_FAIL.smi").open("w", encoding="utf-8")

def props(m):
    return {
        "MW":   Descriptors.MolWt(m),
        "LogP": Descriptors.MolLogP(m),
        "TPSA": rdMolDescriptors.CalcTPSA(m),
        "HBD":  rdMolDescriptors.CalcNumHBD(m),
        "HBA":  rdMolDescriptors.CalcNumHBA(m),
        "ROTB": rdMolDescriptors.CalcNumRotatableBonds(m),
        "AROM": rdMolDescriptors.CalcNumAromaticRings(m),
        "FSP3": Descriptors.FractionCSP3(m),
        "CHARGE": int(sum(a.GetFormalCharge() for a in m.GetAtoms())),
    }

def pass_lenient(p):
    return (200 <= p["MW"] <= 500 and p["TPSA"] <= 90 and p["HBD"] <= 2 and p["HBA"] <= 10
            and 1.5 <= p["LogP"] <= 4.2 and p["ROTB"] <= 10)

def pass_moderate(p):
    return (200 <= p["MW"] <= 450 and p["TPSA"] <= 80 and p["HBD"] <= 2 and p["HBA"] <= 8
            and 1.8 <= p["LogP"] <= 4.0 and p["ROTB"] <= 8
            and 1 <= p["AROM"] <= 4 and 0.15 <= p["FSP3"] <= 0.65
            and p["CHARGE"] >= 0)  # avoid anions

def pass_strict(p):
    return (250 <= p["MW"] <= 430 and p["TPSA"] <= 70 and p["HBD"] <= 1 and p["HBA"] <= 7
            and 2.0 <= p["LogP"] <= 3.5 and p["ROTB"] <= 7
            and 1 <= p["AROM"] <= 3 and 0.20 <= p["FSP3"] <= 0.60
            and p["CHARGE"] >= 0)  # neutral or +1

count = kept_len = kept_mod = kept_str = 0
with open(INPUT_SMI, "r", encoding="utf-8", errors="ignore") as fin:
    for line in fin:
        s = line.strip()
        if not s: continue
        parts = s.split()
        smi = parts[0]
        name = parts[1] if len(parts) > 1 else ""
        m = Chem.MolFromSmiles(smi)
        if m is None:
            f_fail.write(s + "\n"); continue
        p = props(m)

        wrote = False
        if pass_strict(p):
            f_strict.write(s + "\n"); kept_str += 1; wrote = True
        elif pass_moderate(p):
            f_moderate.write(s + "\n"); kept_mod += 1; wrote = True
        elif pass_lenient(p):
            f_lenient.write(s + "\n"); kept_len += 1; wrote = True
        else:
            f_fail.write(s + "\n")
        count += 1
        if count % 1_000_000 == 0:
            print(f"Processed {count:,} | strict {kept_str:,} | moderate {kept_mod:,} | lenient {kept_len:,}")

for fh in (f_strict, f_moderate, f_lenient, f_fail): fh.close()
print(f"[DONE] Processed {count:,}")
print(f"  STRICT:   {kept_str:,}")
print(f"  MODERATE: {kept_mod:,}")
print(f"  LENIENT:  {kept_len:,}")
print("Outputs:")
print(" ", Path(OUT_DIR, 'CNS_STRICT.smi'))
print(" ", Path(OUT_DIR, 'CNS_MODERATE.smi'))
print(" ", Path(OUT_DIR, 'CNS_LENIENT.smi'))
print(" ", Path(OUT_DIR, 'CNS_FAIL.smi'))
