import os, shutil
import pandas as pd

base = r"C:\Users\JDC\Desktop\TIMP2 RESEARCH"
out_dir = os.path.join(base, "All_Docked_Poses")
os.makedirs(out_dir, exist_ok=True)

# List your chunk names
chunks = ["strict_0001", "strict_0028", "strict_0055", "strict_0082", "strict_0109", "strict_0120"]

for chunk in chunks:
    csv_path = os.path.join(base, chunk, f"{chunk}_with_zinc.csv")
    vina_dir = os.path.join(base, chunk, chunk, "vina_out")

    if not os.path.exists(csv_path):
        print(f"⚠️ No CSV for {chunk}")
        continue

    df = pd.read_csv(csv_path)

    # Get name columns safely
    if "ligand" not in df.columns or "zinc_id" not in df.columns:
        print(f"⚠️ Columns missing in {chunk}")
        continue

    for row in df.itertuples():
        lig_file = f"site1_{getattr(row, 'ligand', '')}.pdbqt"
        zinc_id = str(getattr(row, 'zinc_id', '')).strip()
        src = os.path.join(vina_dir, lig_file)
        dst = os.path.join(out_dir, f"{zinc_id}.pdbqt")

        if os.path.exists(src):
            shutil.copy(src, dst)
            print(f"✅ Copied {src} → {zinc_id}.pdbqt")
        else:
            print(f"❌ Missing: {src}")
