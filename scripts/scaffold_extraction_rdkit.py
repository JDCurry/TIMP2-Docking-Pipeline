#!/usr/bin/env python3
"""
Scaffold check for specifc sites
Requires: RDKit (https://www.rdkit.org)
Usage:
    python scaffold_extraction_rdkit.py
Outputs:
    - scaffold_default.smi  (Bemis–Murcko, linkers/side chains removed)
    - scaffold_linkers_kept.smi (Murcko + linkers preserved as core)
    - depiction.png (2D image of ligand + both scaffolds)
"""
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import Draw

smiles = "NC(=O)CCC(=O)O[C@@H](c1ccc(C(F)(F)F)cc1)C(F)(F)F"
mol = Chem.MolFromSmiles(smiles)

# Default BM scaffold (rings only, side chains stripped)
scaf = MurckoScaffold.GetScaffoldForMol(mol)
Chem.MolToSmiles(scaf, isomericSmiles=True)
with open("scaffold_default.smi","w") as f:
    f.write(Chem.MolToSmiles(scaf, isomericSmiles=True)+"\tdefault\n")

# "Linkers preserved" variant: treat linkers as part of core (rings + linkers)
# RDKit 2022.09+ offers GetScaffoldForMol with options, but a common approach
# is to build the Murcko "framework" and then recover linkers via ring-edge paths.
# Here we approximate using the "Murcko generic" then add shortest paths from
# ring atoms to attachment points.

from rdkit.Chem.Scaffolds.MurckoScaffold import MakeScaffoldGeneric
from rdkit.Chem.rdmolops import GetShortestPath

ring_info = mol.GetRingInfo()
ring_atoms = set([a for ring in ring_info.AtomRings() for a in ring])

# Attachment points: atoms in side chains directly bonded to ring atoms
attach_atoms = set()
for a in ring_atoms:
    for nb in mol.GetAtomWithIdx(a).GetNeighbors():
        if nb.GetIdx() not in ring_atoms:
            attach_atoms.add(nb.GetIdx())

# Find linker atoms along shortest paths between ring atoms and first hetero acceptor/donor atoms
hetero = set([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() not in (1,6)])

linker_atoms = set()
for att in attach_atoms:
    # Extend to nearest hetero atom to keep pharmacophore-carrying linker
    targets = list(hetero - ring_atoms)
    if not targets:
        continue
    # Choose nearest hetero atom by topological distance
    d, path_best = 1e9, None
    for t in targets:
        path = GetShortestPath(mol, att, t)
        if path and 0 < len(path) < d:
            d = len(path)
            path_best = path
    if path_best:
        linker_atoms.update(path_best)

core_atoms = sorted(ring_atoms | linker_atoms)
em = Chem.EditableMol(Chem.MolFromSmiles(""))
idx_map = {}
for idx in core_atoms:
    at = mol.GetAtomWithIdx(idx)
    na = Chem.Atom(at.GetAtomicNum())
    na.SetIsAromatic(at.GetIsAromatic())
    new_idx = em.AddAtom(na)
    idx_map[idx] = new_idx

# Add bonds within core
for b in mol.GetBonds():
    a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    if a1 in idx_map and a2 in idx_map:
        em.AddBond(idx_map[a1], idx_map[a2], b.GetBondType())

core = em.GetMol()
Chem.SanitizeMol(core)
core = Chem.RemoveHs(core)
scaf_linkers = Chem.MolFromSmiles(Chem.MolToSmiles(core, isomericSmiles=True))

with open("scaffold_linkers_kept.smi","w") as f:
    f.write(Chem.MolToSmiles(scaf_linkers, isomericSmiles=True)+"\tlinkers_kept\n")

# Draw the ligand and both scaffolds for quick comparison
mols = [mol, scaf, scaf_linkers]
legends = ["Ligand", "Murcko (rings only)", "Murcko+linkers"]
img = Draw.MolsToImage(mols, molsPerRow=3, subImgSize=(300,300), legends=legends)
img.save("depiction.png")
print("Wrote: scaffold_default.smi, scaffold_linkers_kept.smi, depiction.png")
