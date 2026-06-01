#!/usr/bin/env python3
"""
TIMP2 Ligand-Bound MD — Extended from timp2_md_local.py
=========================================================
Runs protein+ligand MD for all 3 Tier 1A ZINC compounds from the paper.
Preserves all original parameters from the working 50 ns run.

Key parameters (unchanged from original):
    padding      = 1.2 nm       ← NOT 80 Å (that caused 762K atom box)
    timestep     = 4.0 fs       (hydrogen mass repartitioning)
    HMR          = 1.5 amu
    forcefield   = AMBER14 + TIP3P-FB
    temperature  = 300 K, 1 bar, 0.15 M NaCl

Ligand parameterization: GAFF2 via openmmforcefields
    No pre-generated MOL2/FRCMOD files needed.
    Works from the PDB files already in PDB_RDKit/.

Hardware target: RTX 6000 (CUDA, with OpenCL fallback)
Expected: ~55-70 ns/day per ligand (~42K atoms)
Wall time: ~35-40 hours per ligand at 100 ns

Dependencies (beyond original):
    pip install openmmforcefields openff-toolkit rdkit

IMPORTANT: Use OpenMM 8.1.2 (not 8.5.1) for CUDA compatibility:
    conda install -c conda-forge openmm=8.1.2

Usage:
    python3 timp2_md_ligand_bound.py --test              # 10 ns, best_binder only (CUDA)
    python3 timp2_md_ligand_bound.py --ns 100            # 100 ns, all 3 ligands (CUDA)
    python3 timp2_md_ligand_bound.py --ns 100 --ligand best_binder
    python3 timp2_md_ligand_bound.py --ns 50 --platform OpenCL
    python3 timp2_md_ligand_bound.py --check             # dep/file check only

Output per ligand (e.g. best_binder_100ns/):
    md_start.pdb              Reference structure
    md_production.dcd         Trajectory
    production.log            Energy/temperature log
    md_analysis.json          Pocket + ligand stability results
    pocket_stability.png      Publication figure (5-panel)
    md_results_latex.tex      Ready-to-paste LaTeX snippet
"""

import os
import sys
import time
import json
import argparse
import math
import warnings
warnings.filterwarnings('ignore')

import numpy as np

# ============================================================
# LIGAND CONFIGURATION
# ============================================================

PROTEIN_PDB = '/home/jd/Desktop/timp2/1br9_fixed.pdb'
LIGAND_DIR  = '/home/jd/Desktop/timp2/PDB_docked'  # Default docked-coordinate directory
NEW_LIGAND_DIR = '/home/jd/Desktop/timp2/best_binder_2'

# Three Tier 1A candidates from the paper
# All ligands now use docked PDB files (Vina coordinates + hydrogens)
LIGANDS = {
    'best_binder': {
        'pdb':    'ZINC000870248439_docked.pdb',
        'zinc':   'ZINC000870248439',
        'label':  'Best Binder (ZINC870248439)',
        'smiles': 'NC(=O)c1ccc(C(=O)N2CCC3(CCSCC3)CC2)cc1',
    },
    'highest_safety': {
        'pdb':    'ZINC000543478179_docked.pdb',
        'zinc':   'ZINC000543478179',
        'label':  'Highest Safety (ZINC543478179)',
        'smiles': 'CCN1CCN(C(=O)N[C@H](C)c2ccc(OCC3CC3)c(F)c2)CC1=O',
    },
    'scaffold_rep': {
        'pdb':    'ZINC000329563636_docked.pdb',
        'zinc':   'ZINC000329563636',
        'label':  'Scaffold Rep (ZINC329563636)',
        'smiles': 'CCN1CCN(C(=O)N[C@H]2CC3(CCCC3)Oc3ccccc32)CC1=O',
    },
    'site1_analog': {
        'pdb':    'ZINC000583127796_docked.pdb',
        'zinc':   'ZINC000583127796',
        'label':  'Site1 Analog (ZINC583127796)',
        'smiles': 'C#Cc1ccc(C(=O)NC[C@@H]2CCCN2C(=O)/C=C(\\C)C2CC2)cn1',
    },
}

# MD parameters — identical to original timp2_md_local.py
MD = {
    'temperature_K':    300,
    'pressure_bar':     1.0,
    'ionic_strength_M': 0.15,
    'timestep_fs':      4.0,       # fs — HMR allows 4 fs safely
    'hmr_amu':          1.5,       # hydrogen mass repartitioning
    'padding_nm':       1.2,       # nm — original used 1.2*nanometers
    'cutoff_nm':        1.0,
    'chunk_ns':         10.0,      # progress report interval
    'save_ps':          50.0,      # trajectory frame interval
    'chk_ns':           5.0,       # checkpoint interval
}

# ============================================================
# ARGUMENT PARSING
# ============================================================

parser = argparse.ArgumentParser(description='TIMP2 Ligand-Bound MD')
parser.add_argument('--ns',       type=float, default=100,
                    help='Production MD length in ns (default: 100)')
parser.add_argument('--test',     action='store_true',
                    help='Quick 10 ns test on best_binder only')
parser.add_argument('--resume',   action='store_true',
                    help='Resume from checkpoint if available')
parser.add_argument('--ligand',   type=str,   default=None,
                    choices=list(LIGANDS.keys()),
                    help='Run a single ligand (default: all 3)')
parser.add_argument('--platform', type=str,   default='CUDA',
                    choices=['CUDA', 'OpenCL', 'CPU'],
                    help='Compute platform (default: CUDA)')
parser.add_argument('--check',    action='store_true',
                    help='Check dependencies and files (including the new ligand path), then exit')
args = parser.parse_args()

if args.test:
    args.ns    = 10.0
    args.ligand = 'best_binder'


# ============================================================
# DEPENDENCY / FILE CHECK
# ============================================================

def check_all():
    ok = True
    print("\nDependency check:")
    for mod, install in [
        ('openmm',            'conda install -c conda-forge openmm'),
        ('openmmforcefields', 'pip install openmmforcefields'),
        ('openff.toolkit',    'pip install openff-toolkit'),
        ('rdkit',             'pip install rdkit'),
        ('pdbfixer',          'conda install -c conda-forge pdbfixer'),
        ('MDAnalysis',        'pip install mdanalysis'),
        ('matplotlib',        'pip install matplotlib'),
    ]:
        try:
            __import__(mod)
            print(f"  ✓ {mod}")
        except ImportError:
            print(f"  ✗ {mod}  →  {install}")
            if mod in ('openmm', 'openmmforcefields', 'openff.toolkit'):
                ok = False  # these are required; others are analysis-only

    print("\nFile check:")
    for label, path in [('Protein PDB', PROTEIN_PDB)] + \
                       [(k, os.path.join(v.get('dir', LIGAND_DIR), v.get('pdb', v.get('pdbqt'))))
                        for k, v in LIGANDS.items()]:
        exists = os.path.exists(path)
        print(f"  {'✓' if exists else '✗'} {label}: {path}")
        if not exists:
            ok = False
    return ok


if args.check:
    sys.exit(0 if check_all() else 1)


# ============================================================
# IMPORTS
# ============================================================

try:
    import openmm
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
except ImportError:
    print("✗ OpenMM not found. Install: conda install -c conda-forge openmm")
    sys.exit(1)

# Plugin path fix (from original script)
import os as _os, sys as _sys
_plugin_dir = _os.environ.get('OPENMM_PLUGIN_DIR') or \
              _os.path.join(_sys.prefix, 'lib', 'plugins')
_os.environ['OPENMM_PLUGIN_DIR'] = _plugin_dir
try:
    openmm.Platform.loadPluginsFromDirectory(_plugin_dir)
except Exception:
    pass

try:
    from openmmforcefields.generators import GAFFTemplateGenerator
except ImportError:
    print("✗ openmmforcefields not found.")
    print("  Install: pip install openmmforcefields")
    sys.exit(1)


def _load_ligand_pose(ligand_cfg, output_dir, ligand_name):
    """Load a docked ligand from PDB or PDBQT; convert PDBQT to PDB if needed."""
    ligand_dir = ligand_cfg.get('dir', LIGAND_DIR)
    ligand_file = ligand_cfg.get('pdb') or ligand_cfg.get('pdbqt')
    ligand_path = os.path.join(ligand_dir, ligand_file)

    if not os.path.exists(ligand_path):
        raise FileNotFoundError(f"Docked ligand not found: {ligand_path}")

    if ligand_path.lower().endswith('.pdbqt'):
        converted_path = os.path.join(output_dir, f'{ligand_name}_ligand.pdb')
        smiles = ligand_cfg.get('smiles')
        atom_records = []
        heavy_atom_records = []
        with open(ligand_path, 'r') as fin, open(converted_path, 'w') as fout:
            for line in fin:
                if line.startswith('REMARK SMILES') and smiles is None:
                    smiles = line.split('REMARK SMILES', 1)[1].strip()
                elif line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = (line[76:78].strip() or atom_name[0]).upper()
                    atom_records.append((atom_name, element, x, y, z))
                    if element != 'H':
                        heavy_atom_records.append((atom_name, element, x, y, z))

        if smiles is None:
            raise ValueError(f'No SMILES string found for {ligand_name} in {ligand_path}')

        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem

            rd_mol = Chem.MolFromSmiles(smiles)
            if rd_mol is None:
                raise ValueError(f'RDKit could not parse SMILES for {ligand_name}')
            rd_mol = Chem.AddHs(rd_mol)

            rd_heavy_atoms = [a for a in rd_mol.GetAtoms() if a.GetSymbol() != 'H']
            if len(rd_heavy_atoms) != len(heavy_atom_records):
                raise ValueError(
                    f'Heavy-atom mismatch for {ligand_name}: '
                    f'PDBQT has {len(heavy_atom_records)} heavy atoms, SMILES has '
                    f'{len(rd_heavy_atoms)} heavy atoms.'
                )

            AllChem.EmbedMolecule(rd_mol, randomSeed=0xC0FFEE)
            conformer = rd_mol.GetConformer()

            heavy_index = 0
            for atom in rd_mol.GetAtoms():
                if atom.GetSymbol() == 'H':
                    continue
                _, _, x, y, z = heavy_atom_records[heavy_index]
                conformer.SetAtomPosition(atom.GetIdx(), (x, y, z))
                heavy_index += 1

            Chem.MolToPDBFile(rd_mol, converted_path)
        except Exception as e:
            raise ValueError(f'Failed to build hydrogenated ligand PDB for {ligand_name}: {e}')

        return PDBFile(converted_path), smiles, ligand_path

    return PDBFile(ligand_path), ligand_cfg.get('smiles'), ligand_path


# ============================================================
# PHASE 1: SYSTEM SETUP (protein + ligand)
# ============================================================

def setup_system(ligand_name, ligand_cfg, output_dir):
    """
    Load protein + docked ligand, parameterize with GAFF2, solvate.
    
    Key fix: GAFFTemplateGenerator requires an openff Molecule object,
    which is created from SMILES (not from a PDB file).
    """
    from openmmforcefields.generators import GAFFTemplateGenerator
    
    solvated_path = os.path.join(output_dir, f'{ligand_name}_solvated.pdb')

    # --- Load protein ---
    print(f"\n  Loading protein: {PROTEIN_PDB}")
    pdb = PDBFile(PROTEIN_PDB)
    print(f"    {pdb.topology.getNumResidues()} residues, "
          f"{pdb.topology.getNumAtoms()} atoms")

    # --- Load docked ligand PDB ---
    lig_pdb, smiles_from_file, ligand_pdb_path = _load_ligand_pose(ligand_cfg, output_dir, ligand_name)
    print(f"  Loading docked ligand: {os.path.basename(ligand_pdb_path)}")
    print(f"    Ligand: {lig_pdb.topology.getNumAtoms()} atoms")

    # --- Set up force field with GAFF2 for ligand ---
    print(f"  Setting up force field (AMBER14 + GAFF2)...")
    ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    
    # Create Molecule from SMILES for GAFFTemplateGenerator
    smiles = ligand_cfg.get('smiles') or smiles_from_file
    if smiles is None:
        raise ValueError(
            f"No SMILES string for {ligand_name}!\n"
            f"Get it from: https://zinc20.docking.org/substances/{ligand_cfg['zinc']}/\n"
            f"Or run: obabel {ligand_pdb_path} -osmi\n"
            f"Then add 'smiles': '...' to the LIGANDS dict."
        )
    
    print(f"  Parameterizing ligand with GAFF2 (from SMILES)...")
    from openff.toolkit import Molecule as OFFMolecule
    off_mol = OFFMolecule.from_smiles(smiles)
    
    # Use the PDB topology directly (it already has correct atoms + hydrogens)
    ligand_topology = lig_pdb.topology
    
    gaff = GAFFTemplateGenerator(molecules=[off_mol], forcefield='gaff-2.11')
    ff.registerTemplateGenerator(gaff.generator)
    print(f"    GAFF2 parameters registered")

    # --- Combine protein + ligand first, then solvate the complex ---
    print(f"  Combining protein + docked ligand...")
    modeller = Modeller(pdb.topology, pdb.positions)
    modeller.add(ligand_topology, lig_pdb.positions)

    print(f"  Solvating complex (1.2 nm padding, 0.15 M NaCl)...")
    modeller.addSolvent(
        ff,
        model='tip3p',
        padding=MD['padding_nm'] * nanometers,
        ionicStrength=MD['ionic_strength_M'] * molar,
    )
    
    n_atoms = modeller.topology.getNumAtoms()
    print(f"    Solvated: {n_atoms:,} atoms")

    PDBFile.writeFile(modeller.topology, modeller.positions,
                      open(solvated_path, 'w'))

    # --- Create system ---
    print(f"  Creating OpenMM system...")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=PME,
        nonbondedCutoff=MD['cutoff_nm'] * nanometers,
        constraints=HBonds,
        hydrogenMass=MD['hmr_amu'] * amu,
    )

    return modeller, system, n_atoms


# ============================================================
# PHASE 2: PLATFORM
# ============================================================

def setup_platform():
    """Configure compute platform — try requested, fallback to fastest available."""
    print(f"\n  Available platforms:")
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        print(f"    {p.getName()}")

    # Try requested platform first
    try:
        platform = Platform.getPlatformByName(args.platform)
        if args.platform == 'CUDA':
            props = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
            print(f"\n  Using CUDA (mixed precision, RTX 6000)")
        elif args.platform == 'OpenCL':
            props = {'Precision': 'mixed'}
            print(f"\n  Using OpenCL (mixed precision)")
        else:
            props = {}
            print(f"\n  Using CPU")
        return platform, props
    except Exception as e:
        print(f"\n  WARNING: {args.platform} unavailable ({e})")

    # Fallback: try OpenCL (usually works on NVIDIA)
    try:
        platform = Platform.getPlatformByName('OpenCL')
        props = {'Precision': 'mixed'}
        print(f"  Falling back to OpenCL (mixed precision)")
        return platform, props
    except Exception:
        pass

    # Final fallback: CPU
    print(f"  Falling back to CPU (slow)")
    platform = Platform.getPlatformByName('CPU')
    props = {}
    return platform, props


def create_simulation_with_fallback(topology, system, integrator, requested_platform, requested_props):
    """
    Try to create Simulation with requested platform.
    On error (e.g., CUDA_ERROR_UNSUPPORTED_PTX_VERSION), auto-fallback to faster available platform.
    """
    platforms_to_try = []

    if requested_platform.getName() == 'CUDA' or requested_platform.getName() == 'OpenCL':
        # If CUDA/OpenCL requested, try both before CPU
        platforms_to_try.append((requested_platform, requested_props, requested_platform.getName()))
        try:
            alt_platform = Platform.getPlatformByName('OpenCL' if requested_platform.getName() == 'CUDA' else 'CUDA')
            alt_props = {'Precision': 'mixed'} if alt_platform.getName() == 'OpenCL' else {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
            platforms_to_try.append((alt_platform, alt_props, alt_platform.getName()))
        except:
            pass
    else:
        platforms_to_try.append((requested_platform, requested_props, requested_platform.getName()))

    # Always add CPU as final fallback
    try:
        cpu = Platform.getPlatformByName('CPU')
        platforms_to_try.append((cpu, {}, 'CPU'))
    except:
        pass

    last_error = None
    for platform, props, name in platforms_to_try:
        try:
            sim = Simulation(topology, system, integrator, platform, props)
            if name != requested_platform.getName():
                print(f"  → Auto-fallback: {name}")
                print(f"     Reason: {last_error}")
            return sim
        except Exception as e:
            last_error = e
            if name == requested_platform.getName() or name != 'CPU':
                print(f"  ✗ {name} failed: {str(e)[:120]}")
            continue

    raise RuntimeError(f"Failed to create Simulation on any available platform: {last_error}")



# ============================================================
# PHASE 3: MINIMIZATION + EQUILIBRATION
# ============================================================

def minimize_and_equilibrate(modeller, system, platform, props, output_dir, ligand_name):
    """Energy minimization, NVT 200 ps, NPT 200 ps."""

    checkpoint_path = os.path.join(output_dir, 'equilibrated.chk')
    start_cif       = os.path.join(output_dir, 'md_start.cif')

    # 4 fs timestep with HMR (same as original)
    integrator = LangevinMiddleIntegrator(
        MD['temperature_K'] * kelvin,
        1 / picosecond,
        MD['timestep_fs'] * 0.001 * picoseconds   # fs → ps
    )

    if args.resume and os.path.exists(checkpoint_path):
        print("\n  Resuming from equilibration checkpoint...")
        simulation = create_simulation_with_fallback(modeller.topology, system, integrator, platform, props)
        simulation.loadCheckpoint(checkpoint_path)
        return simulation

    print("\n[2] Energy minimization...")
    simulation = create_simulation_with_fallback(modeller.topology, system, integrator, platform, props)
    simulation.context.setPositions(modeller.positions)

    e_before = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Initial energy: {e_before}")
    simulation.minimizeEnergy(maxIterations=10000)
    e_after = simulation.context.getState(getEnergy=True).getPotentialEnergy()
    print(f"  Final energy:   {e_after}")

    print("\n[3] NVT equilibration (200 ps, 300 K)...")
    simulation.context.setVelocitiesToTemperature(MD['temperature_K'] * kelvin)
    simulation.reporters.append(
        StateDataReporter(
            os.path.join(output_dir, 'equilibration.log'),
            10000, step=True, temperature=True, potentialEnergy=True, speed=True
        )
    )
    t0 = time.time()
    simulation.step(50000)   # 200 ps at 4 fs
    print(f"  NVT done in {time.time()-t0:.0f}s")

    print("\n[4] NPT equilibration (200 ps, 300 K, 1 bar)...")
    system.addForce(MonteCarloBarostat(1 * bar, MD['temperature_K'] * kelvin))
    simulation.context.reinitialize(preserveState=True)
    t0 = time.time()
    simulation.step(50000)   # 200 ps at 4 fs
    print(f"  NPT done in {time.time()-t0:.0f}s")

    simulation.saveCheckpoint(checkpoint_path)
    state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
    PDBxFile.writeFile(simulation.topology, state.getPositions(),
                       open(start_cif, 'w'))
    print(f"  Saved reference structure: {start_cif}")
    simulation.reporters = []

    return simulation


# ============================================================
# PHASE 4: PRODUCTION MD
# ============================================================

def run_production(simulation, ns, output_dir, ligand_name):
    """Run production MD — same progress format as original."""

    dcd_path = os.path.join(output_dir, 'md_production.dcd')
    log_path = os.path.join(output_dir, 'production.log')
    chk_path = os.path.join(output_dir, 'production.chk')

    dt_fs          = MD['timestep_fs']
    nsteps         = int(ns * 1e6 / dt_fs)
    save_ps        = MD['save_ps'] if ns >= 50 else 10.0
    save_interval  = int(save_ps * 1000 / dt_fs)
    chk_interval   = int(MD['chk_ns'] * 1e6 / dt_fs)
    chunk_steps    = int(MD['chunk_ns'] * 1e6 / dt_fs)

    print(f"\n[5] Production MD ({ns} ns) — {ligand_name}")
    print(f"  Timestep:        {dt_fs} fs (hydrogen mass repartitioning)")
    print(f"  Total steps:     {nsteps:,}")
    print(f"  Save interval:   {save_ps} ps ({save_interval} steps)")
    print(f"  Checkpoint:      every {MD['chk_ns']:.0f} ns")

    resumed = False
    start_step = 0
    if args.resume and os.path.exists(chk_path):
        print(f"\n  Resuming production from checkpoint: {chk_path}")
        simulation.loadCheckpoint(chk_path)
        resumed = True
        start_step = int(getattr(simulation, 'currentStep', 0))
        print(f"  Resumed at step {start_step:,} / {nsteps:,}")

    remaining_steps = max(0, nsteps - start_step)
    n_frames = max(1, math.ceil(remaining_steps / save_interval)) if remaining_steps > 0 else 0
    n_chunks = max(1, math.ceil(remaining_steps / chunk_steps)) if remaining_steps > 0 else 0

    if remaining_steps == 0:
        print("  Nothing left to run; the target trajectory length is already complete.")
        return dcd_path

    simulation.reporters.append(DCDReporter(dcd_path, save_interval, append=resumed and os.path.exists(dcd_path)))
    simulation.reporters.append(
        StateDataReporter(log_path, save_interval,
                          step=True, time=True, temperature=True,
                          potentialEnergy=True, totalEnergy=True,
                          speed=True, remainingTime=True, totalSteps=nsteps,
                          append=resumed and os.path.exists(log_path)))
    simulation.reporters.append(CheckpointReporter(chk_path, chk_interval))

    print(f"\n  Starting production MD...")
    print(f"  {'Progress':>12} {'Time':>8} {'Speed':>14} {'ETA':>10}")
    print(f"  {'-'*48}")

    t_start      = time.time()
    completed_ns = start_step * dt_fs / 1e6

    for i in range(n_chunks):
        batch = min(chunk_steps, remaining_steps - i * chunk_steps)
        if batch <= 0:
            break
        simulation.step(batch)
        completed_ns += MD['chunk_ns']

        elapsed        = time.time() - t_start
        speed_ns_day   = (completed_ns / elapsed) * 86400 if elapsed > 0 else 0
        remaining_ns   = ns - completed_ns
        eta_h          = (remaining_ns / speed_ns_day) * 24 if speed_ns_day > 0 else 0
        pct            = min(100.0, completed_ns / ns * 100)

        print(f"  {pct:>10.1f}%  {completed_ns:>6.1f}ns  "
              f"{speed_ns_day:>10.1f} ns/day  {eta_h:>7.1f}h")

    total_h     = (time.time() - t_start) / 3600
    final_speed = ((ns - start_step * dt_fs / 1e6) / (time.time() - t_start)) * 86400

    print(f"\n  Production MD complete!")
    print(f"  Wall time:   {total_h:.1f} hours")
    print(f"  Performance: {final_speed:.0f} ns/day")
    print(f"  Trajectory:  {dcd_path}")

    return dcd_path


# ============================================================
# PHASE 5: POCKET + LIGAND ANALYSIS
# ============================================================

def analyze_pocket_and_ligand(output_dir, ligand_name, ns, topology=None):
    """
    Extended from original analyze_pocket().
    Adds ligand RMSD and pocket-ligand distance tracking
    on top of the original triad geometry analysis.
    """
    try:
        import MDAnalysis as mda
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print("  ⚠ MDAnalysis/matplotlib not installed — skipping analysis")
        print("    Install: pip install mdanalysis matplotlib")
        return None

    dcd_path = os.path.join(output_dir, 'md_production.dcd')

    topology_source = topology
    if topology_source is None:
        cif_path = os.path.join(output_dir, 'md_start.cif')
        pdb_path = os.path.join(output_dir, 'md_start.pdb')
        topology_source = cif_path if os.path.exists(cif_path) else pdb_path

    print(f"\n[6] Analyzing pocket + ligand stability...")

    u = mda.Universe(topology_source, dcd_path)
    n_frames = len(u.trajectory)
    print(f"  Loaded: {n_frames} frames, {ns:.0f} ns")

    # ---- Pocket triad (Val6, Leu100, Phe103 Cα) ----
    # These are the residues from the paper
    try:
        triad_sel = u.select_atoms('protein and name CA and (resid 6 or resid 100 or resid 103)')
        if len(triad_sel) != 3:
            raise ValueError(f"Expected 3 triad atoms, got {len(triad_sel)}")
    except Exception as e:
        print(f"  ⚠ Triad selection failed: {e}")
        print("    Check residue numbers match 1BR9 numbering.")
        triad_sel = None

    # ---- Protein backbone for RMSD ----
    protein_bb = u.select_atoms('protein and backbone')
    pocket_res  = u.select_atoms(
        'protein and (resid 6 or resid 8 or resid 96 or resid 97 or '
        'resid 99 or resid 100 or resid 103 or resid 104)'
    )
    # ---- Ligand ----
    ligand_sel = u.select_atoms(
        'not protein and not resname HOH WAT TIP3 TIP3P SOL NA CL K'
    )
    has_ligand = len(ligand_sel) > 0
    if has_ligand:
        print(f"  Ligand atoms found: {len(ligand_sel)}")
    else:
        print(f"  ⚠ No ligand atoms found (resname LIG/UNL/MOL) — "
              f"check residue name in md_start.pdb")

    # ---- Reference frame 0 for RMSD ----
    ref_pos_bb     = protein_bb.positions.copy()
    ref_pos_pocket = pocket_res.positions.copy()
    ref_pos_lig    = ligand_sel.positions.copy() if has_ligand else None

    times          = []
    d_6_100        = []
    d_6_103        = []
    d_100_103      = []
    areas          = []
    pocket_rmsds   = []
    protein_rmsds  = []
    ligand_rmsds   = []
    lig_pocket_dist= []   # ligand COM to pocket COM

    for ts in u.trajectory:
        t_ns = ts.time / 1000.0
        times.append(t_ns)

        # Protein RMSD
        protein_rmsds.append(
            np.sqrt(np.mean(np.sum((protein_bb.positions - ref_pos_bb)**2, axis=1)))
        )
        # Pocket RMSD
        pocket_rmsds.append(
            np.sqrt(np.mean(np.sum((pocket_res.positions - ref_pos_pocket)**2, axis=1)))
        )

        # Triad distances + area
        if triad_sel is not None and len(triad_sel) == 3:
            p1, p2, p3 = triad_sel.positions
            d12 = np.linalg.norm(p1 - p2)
            d13 = np.linalg.norm(p1 - p3)
            d23 = np.linalg.norm(p2 - p3)
            d_6_100.append(d12)
            d_6_103.append(d13)
            d_100_103.append(d23)
            # Heron's formula for triangle area
            s = (d12 + d13 + d23) / 2
            arg = max(0.0, s * (s-d12) * (s-d13) * (s-d23))
            areas.append(np.sqrt(arg))

        # Ligand RMSD + pocket-ligand distance
        if has_ligand and ref_pos_lig is not None:
            lig_com     = ligand_sel.center_of_mass()
            pocket_com  = pocket_res.center_of_mass()
            lig_pocket_dist.append(np.linalg.norm(lig_com - pocket_com))
            ligand_rmsds.append(
                np.sqrt(np.mean(np.sum((ligand_sel.positions - ref_pos_lig)**2, axis=1)))
            )

    times = np.array(times)

    # ---- Crystal reference values (from 1BR9 structure) ----
    # These match the values used in the original paper analysis
    crystal_d_6_100   = 11.98  # Å (Val6–Leu100)
    crystal_d_6_103   = 8.22   # Å (Val6–Phe103)
    crystal_d_100_103 = 7.61   # Å (Leu100–Phe103)
    crystal_area      = None
    if crystal_d_6_100 and crystal_d_6_103 and crystal_d_100_103:
        s = (crystal_d_6_100 + crystal_d_6_103 + crystal_d_100_103) / 2
        arg = max(0.0, s * (s-crystal_d_6_100) * (s-crystal_d_6_103) * (s-crystal_d_100_103))
        crystal_area = np.sqrt(arg)

    # ---- Figure ----
    n_rows = 3 if has_ligand else 2
    fig    = plt.figure(figsize=(16, 5 * n_rows))
    gs     = GridSpec(n_rows, 2, figure=fig, hspace=0.45, wspace=0.35)

    # Panel A: Triad distances
    if d_6_100:
        ax = fig.add_subplot(gs[0, 0])
        win = max(1, len(times) // 50)
        for d, label, color, cryst in [
            (d_6_100,   'Val6–Leu100',  '#1f77b4', crystal_d_6_100),
            (d_6_103,   'Val6–Phe103',  '#ff7f0e', crystal_d_6_103),
            (d_100_103, 'Leu100–Phe103','#2ca02c', crystal_d_100_103),
        ]:
            ax.plot(times, d, alpha=0.25, color=color, lw=0.8)
            smooth = np.convolve(d, np.ones(win)/win, mode='same')
            ax.plot(times, smooth, color=color, lw=2, label=label)
            ax.axhline(cryst, color=color, ls='--', alpha=0.6, lw=1)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('(A) Hydrophobic Triad Distances', fontweight='bold')
        ax.legend(fontsize=9)

    # Panel B: Pocket area
    if areas:
        ax = fig.add_subplot(gs[0, 1])
        ax.plot(times, areas, alpha=0.3, color='purple', lw=0.8)
        win = max(1, len(times) // 50)
        smooth = np.convolve(areas, np.ones(win)/win, mode='same')
        ax.plot(times, smooth, color='purple', lw=2)
        if crystal_area:
            ax.axhline(crystal_area, color='red', ls='--', lw=1.5,
                       label=f'Crystal ({crystal_area:.1f} Å²)')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Area (Å²)')
        ax.set_title('(B) Pocket Triangle Area', fontweight='bold')
        ax.legend(fontsize=9)

    # Panel C: RMSD
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(times, protein_rmsds, lw=1.5, color='steelblue', label='Protein backbone')
    ax.plot(times, pocket_rmsds,  lw=1.5, color='darkorange', label='Pocket residues')
    ax.set_xlabel('Time (ns)')
    ax.set_ylabel('RMSD (Å)')
    ax.set_title('(C) RMSD', fontweight='bold')
    ax.legend(fontsize=9)

    # Panel D: Triad distance distributions
    if d_6_100:
        ax = fig.add_subplot(gs[1, 1])
        parts = ax.violinplot(
            [d_6_100, d_6_103, d_100_103],
            positions=[1, 2, 3], showmeans=True, showmedians=True
        )
        for pc, color in zip(parts['bodies'],
                              ['#1f77b4', '#ff7f0e', '#2ca02c']):
            pc.set_facecolor(color); pc.set_alpha(0.5)
        ax.scatter([1, 2, 3],
                   [crystal_d_6_100, crystal_d_6_103, crystal_d_100_103],
                   color='red', s=80, zorder=5, marker='D',
                   label='Crystal structure')
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(['Val6–Leu100', 'Val6–Phe103', 'Leu100–Phe103'],
                            fontsize=10)
        ax.set_ylabel('Distance (Å)')
        ax.set_title('(D) Distance Distributions', fontweight='bold')
        ax.legend(fontsize=9)

    # Panels E–F: Ligand-specific (new)
    if has_ligand and ligand_rmsds:
        ax = fig.add_subplot(gs[2, 0])
        ax.plot(times, ligand_rmsds, color='crimson', lw=1.5)
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('RMSD (Å)')
        ax.set_title('(E) Ligand RMSD', fontweight='bold')

        ax = fig.add_subplot(gs[2, 1])
        ax.plot(times, lig_pocket_dist, color='darkgreen', lw=1.5)
        ax.axhline(np.mean(lig_pocket_dist), color='darkgreen', ls='--',
                   lw=1, label=f'Mean {np.mean(lig_pocket_dist):.1f} Å')
        ax.set_xlabel('Time (ns)')
        ax.set_ylabel('Distance (Å)')
        ax.set_title('(F) Ligand COM – Pocket COM Distance', fontweight='bold')
        ax.legend(fontsize=9)

    lig_label = LIGANDS[ligand_name]['label']
    fig.suptitle(
        f'TIMP2 Cryptic Pocket + Ligand Stability\n'
        f'{ligand_name}: {lig_label}\n'
        f'{times[-1]:.0f} ns MD (AMBER14/TIP3P-FB, 300 K, 1 bar, 0.15 M NaCl)',
        fontsize=13, fontweight='bold', y=0.99
    )

    png_path = os.path.join(output_dir, 'pocket_stability.png')
    pdf_path = os.path.join(output_dir, 'pocket_stability.pdf')
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {png_path}")

    # ---- Verdict ----
    if areas:
        area_preserved = abs(np.mean(areas) - crystal_area) / crystal_area
        cv             = np.std(areas) / np.mean(areas)
        if area_preserved < 0.20 and cv < 0.30:
            verdict = "STABLE"
        elif area_preserved < 0.40 or cv < 0.40:
            verdict = "PARTIALLY STABLE / DYNAMIC"
        else:
            verdict = "UNSTABLE"

        print(f"\n  {'='*55}")
        print(f"  POCKET STABILITY: {verdict}")
        print(f"  Area preservation: {(1-area_preserved)*100:.0f}%")
        print(f"  Coefficient of variation: {cv*100:.0f}%")
    else:
        verdict = "UNKNOWN (triad selection failed)"

    if has_ligand and lig_pocket_dist:
        lig_stable = np.mean(lig_pocket_dist) < 8.0 and np.std(lig_pocket_dist) < 3.0
        print(f"  LIGAND STABILITY: "
              f"{'BOUND' if lig_stable else 'DRIFTING/UNBOUND'}")
        print(f"  Pocket-ligand distance: "
              f"{np.mean(lig_pocket_dist):.1f} ± {np.std(lig_pocket_dist):.1f} Å")
        print(f"  {'='*55}")

    # ---- Save JSON ----
    results = {
        'ligand': ligand_name,
        'simulation': {
            'length_ns':          float(times[-1]),
            'n_frames':           len(times),
            'force_field':        'AMBER14-all + TIP3P-FB + GAFF2',
            'temperature_K':      300,
            'pressure_bar':       1.0,
            'ionic_strength_M':   0.15,
            'timestep_fs':        MD['timestep_fs'],
            'hmr':                True,
        },
        'pocket': {
            'verdict': verdict,
            'triad_distances': {
                'Val6_Leu100':  {'mean': float(np.mean(d_6_100))  if d_6_100  else None,
                                 'std':  float(np.std(d_6_100))   if d_6_100  else None,
                                 'crystal': crystal_d_6_100},
                'Val6_Phe103':  {'mean': float(np.mean(d_6_103))  if d_6_103  else None,
                                 'std':  float(np.std(d_6_103))   if d_6_103  else None,
                                 'crystal': crystal_d_6_103},
                'Leu100_Phe103':{'mean': float(np.mean(d_100_103))if d_100_103 else None,
                                 'std':  float(np.std(d_100_103)) if d_100_103 else None,
                                 'crystal': crystal_d_100_103},
            },
        },
        'rmsd': {
            'pocket_mean':   float(np.mean(pocket_rmsds)),
            'protein_mean':  float(np.mean(protein_rmsds)),
        },
    }

    if has_ligand and ligand_rmsds:
        results['ligand_rmsd'] = [float(value) for value in ligand_rmsds]
        results['pocket_ligand_distance'] = [float(value) for value in lig_pocket_dist]
        results['ligand_dynamics'] = {
            'rmsd_mean_A':          float(np.mean(ligand_rmsds)),
            'rmsd_std_A':           float(np.std(ligand_rmsds)),
            'pocket_distance_mean': float(np.mean(lig_pocket_dist)),
            'pocket_distance_std':  float(np.std(lig_pocket_dist)),
            'bound_throughout':     bool(np.mean(lig_pocket_dist) < 8.0),
        }

    json_path = os.path.join(output_dir, 'md_analysis.json')
    with open(json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\n  Analysis saved: {json_path}")

    return results


# ============================================================
# MAIN — run all 3 ligands sequentially
# ============================================================

if __name__ == '__main__':

    ligands_to_run = ([args.ligand] if args.ligand
                      else list(LIGANDS.keys()))

    print("=" * 65)
    print("TIMP2 Ligand-Bound MD Simulation")
    print(f"Production length: {args.ns} ns per ligand")
    print(f"Platform: {args.platform}")
    print(f"Ligands: {', '.join(ligands_to_run)}")
    print("=" * 65)

    # One-time platform setup
    platform, props = setup_platform()

    summary = []

    for ligand_name in ligands_to_run:
        cfg        = LIGANDS[ligand_name]
        output_dir = os.path.join(
            '/home/jd/Desktop/timp2',
            f'{ligand_name}_{args.ns:.0f}ns'
        )
        os.makedirs(output_dir, exist_ok=True)

        lig_pdb = os.path.join(cfg.get('dir', LIGAND_DIR), cfg.get('pdb', cfg.get('pdbqt')))
        if not os.path.exists(lig_pdb):
            print(f"\n✗ Ligand PDB not found: {lig_pdb}")
            summary.append({'ligand': ligand_name, 'status': 'SKIPPED — PDB not found'})
            continue

        print(f"\n{'='*65}")
        print(f"LIGAND: {ligand_name} — {cfg['label']}")
        print(f"Output: {output_dir}")
        print(f"{'='*65}")

        try:
            print(f"\n[1] Preparing system...")
            modeller, system, n_atoms = setup_system(ligand_name, cfg, output_dir)

            simulation = minimize_and_equilibrate(
                modeller, system, platform, props, output_dir, ligand_name
            )

            t_run_start = time.time()
            run_production(simulation, args.ns, output_dir, ligand_name)
            wall_h = (time.time() - t_run_start) / 3600

            results = analyze_pocket_and_ligand(output_dir, ligand_name, args.ns, simulation.topology)

            summary.append({
                'ligand':       ligand_name,
                'status':       'COMPLETE',
                'wall_hours':   round(wall_h, 1),
                'ns_per_day':   round(args.ns / wall_h * 24, 0) if wall_h > 0 else 0,
                'n_atoms':      n_atoms,
                'verdict':      results.get('pocket', {}).get('verdict', '—')
                                if results else '—',
            })

        except Exception as e:
            import traceback
            traceback.print_exc()
            summary.append({'ligand': ligand_name, 'status': f'ERROR: {e}'})
            continue

    # Final summary
    print(f"\n{'='*65}")
    print(f"SUMMARY")
    print(f"{'='*65}")
    for s in summary:
        if s['status'] == 'COMPLETE':
            print(f"  ✓ {s['ligand']:<20} {s['ns_per_day']:.0f} ns/day  "
                  f"{s['wall_hours']:.1f}h  {s['verdict']}")
        else:
            print(f"  ✗ {s['ligand']:<20} {s['status']}")
    print(f"{'='*65}")
