#!/usr/bin/env python3
"""
TIMP2 Ligand-Bound MD (100 ns) - Optimized for Pre-solvated Complexes
Uses OpenCL/CPU for GPU compatibility (RTX 6000 support)
"""

import os
import sys
import time
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

import numpy as np

# Parse arguments
import argparse
parser = argparse.ArgumentParser(description='TIMP2 Ligand-Bound MD (100 ns)')
parser.add_argument('--ns', type=float, default=100, help='Simulation length (ns)')
parser.add_argument('--test', action='store_true', help='Quick test (5 ns)')
parser.add_argument('--platform', type=str, default='OpenCL', choices=['OpenCL','CPU','CUDA'], help='Platform')
args = parser.parse_args()

if args.test:
    args.ns = 5

LIGANDS = [
    'ZINC000870248439_1',
    'ZINC000543478179_1', 
    'ZINC000329563636_1',
]

print("=" * 70)
print(f"TIMP2 Ligand-Bound MD ({args.ns} ns)")
print(f"Platform: {args.platform}")
print("=" * 70)

# Import OpenMM
try:
    import openmm
    from openmm.app import *
    from openmm import *
    from openmm.unit import *
    print(f"  ✓ OpenMM installed")
except ImportError as e:
    print(f"  ✗ ERROR: {e}")
    sys.exit(1)

# Platform
try:
    platform = Platform.getPlatformByName(args.platform)
    if args.platform == 'OpenCL':
        props = {'Precision': 'mixed'}
    elif args.platform == 'CUDA':
        props = {'CudaPrecision': 'mixed', 'DeviceIndex': '0'}
    else:
        props = {}
    print(f"  ✓ Platform: {platform.getName()}")
except Exception as e:
    print(f"  Falling back to CPU: {e}")
    platform = Platform.getPlatformByName('CPU')
    props = {}

def run_md_simulation(ligand_name, ns):
    """Run MD for one ligand."""
    
    solvated_pdb = os.path.join('PDB_RDKit', f'{ligand_name}_solvated.pdb')
    
    if not os.path.exists(solvated_pdb):
        print(f"  ✗ Solvated PDB not found: {solvated_pdb}")
        return False
    
    print(f"\n{'='*70}")
    print(f"MD: {ligand_name}")
    print(f"{'='*70}")
    
    # Load solvated system
    print(f"  → Loading solvated structure...")
    pdb = PDBFile(solvated_pdb)
    n_atoms = pdb.topology.getNumAtoms()
    print(f"    Atoms: {n_atoms}")
    
    # Create forcefield and system
    print(f"  → Creating system...")
    ff = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    system = ff.createSystem(pdb.topology,
                            nonbondedMethod=PME,
                            nonbondedCutoff=1.0*nanometer,
                            constraints=HBonds,
                            hydrogenMass=1.5*amu)
    
    # Setup integrator and simulation
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 4*femtoseconds)
    system.addForce(MonteCarloBarostat(1*bar, 300*kelvin, 25))
    
    simulation = Simulation(pdb.topology, system, integrator, platform, props)
    simulation.context.setPositions(pdb.positions)
    
    # Minimize
    print(f"  → Energy minimization...")
    simulation.minimizeEnergy(maxIterations=10000)
    
    # Equilibrate
    print(f"  → Equilibration (200 ps)...")
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    simulation.step(50000)  # 200 ps at 4fs
    
    # Production MD
    print(f"  → Production MD ({ns} ns)...")
    dcd_file = os.path.join('PDB_RDKit', f'lig_{ligand_name}_{int(ns)}ns.dcd')
    log_file = os.path.join('PDB_RDKit', f'lig_{ligand_name}_{int(ns)}ns.log')
    
    dt_fs = 4.0
    nsteps = int(ns * 1e6 / dt_fs)
    save_interval = int(50 * 1000 / dt_fs)  # 50ps
    
    simulation.reporters.append(DCDReporter(dcd_file, save_interval))
    simulation.reporters.append(StateDataReporter(log_file, save_interval,
                                                   step=True, time=True, 
                                                   temperature=True,
                                                   potentialEnergy=True,
                                                   totalEnergy=True,
                                                   speed=True))
    
    # Run in chunks
    chunk_ns = 10.0
    chunk_steps = int(chunk_ns * 1e6 / dt_fs)
    n_chunks = int(ns / chunk_ns)
    
    print(f"    {'Progress':>10} {'Time':>8} {'Speed':>12} {'ETA':>8}")
    print(f"    {'-'*45}")
    
    t0 = time.time()
    for i in range(n_chunks):
        completed_ns = i * chunk_ns
        simulation.step(chunk_steps)
        
        elapsed = time.time() - t0
        speed = (completed_ns / elapsed) * 86400 if elapsed > 0 else 0
        remaining = ns - completed_ns
        eta = (remaining / speed) * 24 if speed > 0 else 0
        
        pct = min(100, completed_ns / ns * 100)
        print(f"    {pct:>9.0f}%  {completed_ns:>6.1f}ns  {speed:>8.0f} ns/day  {eta:>6.1f}h")
    
    elapsed = time.time() - t0
    print(f"\n  ✓ Complete in {elapsed/3600:.1f} hours")
    print(f"    Trajectory: {os.path.basename(dcd_file)}")
    
    return True

if __name__ == '__main__':
    t_start = time.time()
    
    print(f"\n[START] {datetime.now().strftime('%H:%M:%S')}")
    
    success_count = 0
    for ligand in LIGANDS:
        try:
            if run_md_simulation(ligand, args.ns):
                success_count += 1
        except Exception as e:
            print(f"  ✗ Error: {e}")
            import traceback
            traceback.print_exc()
    
    elapsed_total = time.time() - t_start
    print(f"\n{'='*70}")
    print(f"COMPLETE: {success_count}/{len(LIGANDS)} simulations finished")
    print(f"Total time: {elapsed_total/3600:.1f} hours")
    print(f"End: {datetime.now().strftime('%H:%M:%S')}")
    print(f"{'='*70}")
