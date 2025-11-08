import os

cmd.delete("all")

# === Load TIMP2 receptor and prep ===
cmd.fetch("1BR9", "TIMP2")
cmd.remove("solvent")
cmd.remove("TIMP2 and not chain A")
cmd.hide("everything")
cmd.show("cartoon", "TIMP2")
cmd.color("gray80", "TIMP2")

# === Define active site pocket ===
cmd.pseudoatom("site2_center", pos=[23.255162,27.655522,14.632348], color="salmon")
cmd.select("site2_vicinity", "byres (TIMP2 within 5 of site2_center)")
cmd.show("surface", "site2_vicinity")
cmd.color("palegreen", "site2_vicinity")
cmd.set("transparency", 0.55)

# === Rendering parameters ===
cmd.set("ray_trace_mode", 1)
cmd.set("ray_opaque_background", 0)
cmd.set("antialias", 2)
cmd.set("ray_shadows", 0)

output_dir = r""
os.makedirs(output_dir, exist_ok=True)

# === Ligand list (15) ===
ligs = [
    "ZINC000329563636_1.pdb",
    "ZINC000485614589_1.pdb",
    "ZINC000534546662_1.pdb",
    "ZINC000571415585_1.pdb",
    "ZINC000583053301_1.pdb",
    "ZINC000661273949_1.pdb",
    "ZINC000781056905_1.pdb",
    "ZINC000934141726_1.pdb",
    "ZINC001016817530_1.pdb",
    "ZINC001228612316_1.pdb",
    "ZINC001547745984_1.pdb",
    "ZINC001588320882_1.pdb",
    "ZINC001621861216_1.pdb",
    "ZINC001647470722_1.pdb",
    "ZINC001746475657_1.pdb"
]

# === Render each ligand ===
for lig in ligs:
    name = os.path.splitext(lig)[0]
    full_path = os.path.join(r"C:\", lig)
    print(f"Rendering {name} ...")
    if not os.path.exists(full_path):
        print(f" Missing file: {full_path}")
        continue
    cmd.delete("ligand_temp")
    cmd.load(full_path, "ligand_temp")
    cmd.hide("everything")
    cmd.show("cartoon", "TIMP2")
    cmd.show("surface", "site2_vicinity")
    cmd.color("palegreen", "site2_vicinity")
    cmd.show("sticks", "ligand_temp")
    cmd.util.cbag("ligand_temp")
    cmd.zoom("ligand_temp", 6)
    cmd.ray(2400, 1800)
    out_path = os.path.join(output_dir, f"{name}_site2.png")
    cmd.png(out_path, dpi=300)
    print(f" Saved {out_path}")

print(" All 15 ligand renderings complete.")
