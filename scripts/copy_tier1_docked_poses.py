import os
import shutil
import pandas as pd

# === CONFIG ===
base = r"C:\Users\JDC\Desktop\TIMP2 RESEARCH"
output_dir = os.path.join(base, "tier1_docked_poses")
os.makedirs(output_dir, exist_ok=True)

tier1_csv = os.path.join(base, "File_S3_Tier1_Compounds_SMILES.csv")
tier1 = pd.read_csv(tier1_csv)
allowed_ids = set(tier1["ZINC ID"].astype(str).str.strip())

chunks = [
    "strict_0001",
    "strict_0028",
    "strict_0055",
    "strict_0082",
    "strict_0109",
    "strict_0120"
]

missing_records = []

def find_pdbqt(base_dir, ligand_name):
    """Search recursively for .pdbqt file with matching ligand base name."""
    ligand_core = ligand_name.lower().replace(".pdbqt", "").replace("site1_", "").replace("site2_", "")
    for root, _, files in os.walk(base_dir):
        for f in files:
            if f.lower().endswith(".pdbqt"):
                fname_core = f.lower().replace(".pdbqt", "")
                # Match even if extra prefix/suffix like 'site1_' or '_out'
                if ligand_core in fname_core or fname_core in ligand_core:
                    return os.path.join(root, f)
    return None

# === MAIN LOOP ===
for chunk in chunks:
    chunk_dir = os.path.join(base, chunk)
    csv_path = next((os.path.join(chunk_dir, f)
                     for f in os.listdir(chunk_dir)
                     if f.endswith("_with_zinc.csv")), None)

    if not csv_path:
        print(f"⚠️ Missing mapping CSV for {chunk}")
        continue

    df = pd.read_csv(csv_path)
    df.columns = [c.strip().lower() for c in df.columns]

    ligand_col = "ligand" if "ligand" in df.columns else "id" if "id" in df.columns else None
    zinc_col = "zinc_id" if "zinc_id" in df.columns else "zinc id"

    if not ligand_col:
        print(f"⚠️ No ligand/id column found in {chunk}: {df.columns}")
        continue

    copied = set()
    for _, row in df.iterrows():
        zinc_id = str(row[zinc_col]).strip()
        if zinc_id not in allowed_ids or zinc_id in copied:
            continue

        ligand_name = str(row[ligand_col]).strip()
        expected_file = f"{ligand_name}.pdbqt"

        # First, try directly in vina_out or root
        vina_folder = os.path.join(chunk_dir, "vina_out")
        src_path = os.path.join(vina_folder, expected_file)
        if not os.path.exists(src_path):
            # recursive fallback search
            src_path = find_pdbqt(chunk_dir, ligand_name)

        if src_path and os.path.exists(src_path):
            dest_path = os.path.join(output_dir, f"{zinc_id}_1.pdbqt")
            shutil.copy(src_path, dest_path)
            copied.add(zinc_id)
            print(f"✅ Copied {zinc_id} from {chunk} ({os.path.basename(src_path)})")
        else:
            missing_records.append({"chunk": chunk, "ligand": ligand_name, "zinc_id": zinc_id})
            print(f"❌ Missing file: {ligand_name}.pdbqt in {chunk}")

print("🎉 Done. All Tier-1 docked poses collected.")

# write missing list
if missing_records:
    missing_csv = os.path.join(output_dir, "missing_tier1_files.csv")
    pd.DataFrame(missing_records).to_csv(missing_csv, index=False)
    print(f"📄 Missing files logged to: {missing_csv}")
