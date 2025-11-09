import pandas as pd
from pathlib import Path

# === CONFIGURATION ===
base_dir = Path(r"C:\Users\JDC\Desktop\sites")
out_site1 = base_dir / "site1_scores_master.csv"
out_site2 = base_dir / "site2_scores_master.csv"
out_delta = base_dir / "delta_scores_master.csv"
out_avg = base_dir / "delta_scores_summary.csv"

# === GATHER ALL STRICT FILES ===
csv_files = sorted(base_dir.glob("strict_*_with_zinc.csv"))
print(f"Found {len(csv_files)} chunk files.")

site1_list, site2_list, delta_list = [], [], []

for f in csv_files:
    try:
        df = pd.read_csv(f)

        # check for essential columns
        if not {"site", "score_kcal_per_mol", "ligand", "zinc_id"}.issubset(df.columns):
            print(f"Skipping {f.name} (missing expected columns)")
            continue

        # separate by site
        s1 = (
            df[df["site"] == "site1"]
            .loc[:, ["ligand", "zinc_id", "score_kcal_per_mol"]]
            .assign(chunk=f.stem)
        )
        s2 = (
            df[df["site"] == "site2"]
            .loc[:, ["ligand", "zinc_id", "score_kcal_per_mol"]]
            .assign(chunk=f.stem)
        )

        site1_list.append(s1)
        site2_list.append(s2)

        # compute ΔScore per file
        if not s1.empty and not s2.empty:
            merged = pd.merge(
                s1,
                s2,
                on=["ligand", "zinc_id", "chunk"],
                suffixes=("_site1", "_site2"),
                how="inner",
            )
            merged["delta_score"] = (
                merged["score_kcal_per_mol_site2"] - merged["score_kcal_per_mol_site1"]
            )
            delta_list.append(merged)

        print(f"Processed {f.name}: {len(s1)} site1, {len(s2)} site2 entries")

    except Exception as e:
        print(f"Error reading {f.name}: {e}")

# === CONCATENATE AND SAVE ===
def safe_concat(dfs):
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

site1_master = safe_concat(site1_list)
site2_master = safe_concat(site2_list)
delta_master = safe_concat(delta_list)

site1_master.to_csv(out_site1, index=False)
site2_master.to_csv(out_site2, index=False)
delta_master.to_csv(out_delta, index=False)

# === AVERAGE ΔSCORE SUMMARY ===
if not delta_master.empty:
    avg_summary = (
        delta_master.groupby(["ligand", "zinc_id"], as_index=False)
        .agg(
            avg_site1=("score_kcal_per_mol_site1", "mean"),
            avg_site2=("score_kcal_per_mol_site2", "mean"),
            avg_delta=("delta_score", "mean"),
        )
        .sort_values("avg_delta")
    )
    avg_summary.to_csv(out_avg, index=False)
    print(f"✅ averaged ΔScore summary saved to: {out_avg} ({len(avg_summary)} ligands)")

print(f"\n✅ site1 master: {out_site1}")
print(f"✅ site2 master: {out_site2}")
print(f"✅ delta master: {out_delta}")
