from pathlib import Path
import pandas as pd
import re

# --- Config ---
chunk_dir = Path(r"D:/timp2/work_strict/strict_0028")
vina_dir = chunk_dir / "vina_out"  # or "poses_keep" if you only want the kept poses

rows = []
for p in vina_dir.glob("*.pdbqt"):
    smiles = None
    score = None
    try:
        with open(p, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith("REMARK SMILES") and smiles is None:
                    # grab everything after 'REMARK SMILES'
                    smiles = line.split(maxsplit=2)[2].strip()
                if "REMARK VINA RESULT" in line and score is None:
                    m = re.search(r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)", line)
                    if m:
                        score = float(m.group(1))
            if smiles is not None and score is not None:
                rows.append({"id": p.stem, "docking_score": score, "smiles": smiles})
    except Exception as e:
        print(f"Skipping {p.name} due to error: {e}")

if not rows:
    raise RuntimeError(f"No ligands parsed from {vina_dir}")

hits = pd.DataFrame(rows)
hits.to_csv(chunk_dir / "strict_0028_hits.csv", index=False)
print(f"Recovered {len(hits)} ligands with SMILES and scores from {vina_dir}.")