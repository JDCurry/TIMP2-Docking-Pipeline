# TIMP2 Docking and ADMET Pipeline
*A computational framework for identifying novel allosteric modulators of TIMP2 and exploring mechanisms of neuroplasticity.*

---

## Overview
This repository contains all custom Python scripts, documentation, and supporting materials for the study:

> **Curry, J.D. (2025).**  
> *Large-Scale Virtual Screening Identifies Novel Allosteric Modulators of TIMP2:  
> A Computational Approach to Enhancing Neuroplasticity.*  
> (Preprint forthcoming on OSF)

TIMP2 (Tissue Inhibitor of Metalloproteinase-2) regulates hippocampal plasticity and cognitive function through MMP-independent neuronal signaling. Building on Ferreira *et al.* (2023, *eNeuro*):contentReference[oaicite:0]{index=0}, this work identifies a cryptic, chemically promiscuous allosteric pocket in TIMP2 and demonstrates that it accommodates both halogenated and non-halogenated scaffolds with high CNS drug-likeness.

---

## Biological Context
TIMP2 enhances neuroplasticity by stabilizing neuronal ECM dynamics and engaging α3β1 integrin–Rap1–ERK pathways independently of metalloproteinase inhibition. Small-molecule allosteric modulation offers a potential therapeutic route to replicate these effects with greater bioavailability and blood-brain barrier penetration.

Reference:  
Ferreira T.A. *et al.* (2023). *Neuronal TIMP2 Regulates Hippocampal Neurogenesis and Synaptic Plasticity.* *eNeuro*, 10(4): ENEURO.0031-23.2023.

---

## ⚙️ Pipeline Summary

### Primary Workflow
| Step | Script | Function |
|------|---------|-----------|
| 1️⃣ | `dock_adaptive_fixed.py` | Adaptive AutoDock Vina screening with CPU load benchmarking. |
| 2️⃣ | `rebuild_hits.py` | Parses `.pdbqt` files to reconstruct SMILES and docking scores. |
| 3️⃣ | `timp2_triage.py` | Triage CNS/peripheral hits; filters based on energy, site, and diversity. |
| 4️⃣ | `timp2_analysis.py` | Integrates ADMET data, calculates CNS MPO, clusters scaffolds. |
| 5️⃣ | `add_zinc_ids_prioritized.py` | Re-maps prioritized compounds to their original ZINC IDs. |
| 6️⃣ | `tier_1_filter.py` | Applies strict ADMET-AI toxicity flags (Tier 1 = clean). |
| 7️⃣ | `diversity_selection_script.py` | Balances chemical diversity across chunks. |
| 8️⃣ | `singleton.py` | Isolates unique scaffolds (singleton vs non-singleton). |
| 9️⃣ | `scaffold_extraction_rdkit.py` | Generates Murcko scaffold images for visual QC. |

---

## Computational Setup

**Hardware Used:**
- Dual Xeon workstation (24 cores, 3.2 GHz, 128 GB RAM)
- NVMe SSD, Windows 10 64-bit
- AutoDock Vina 1.2.5, RDKit 2023.03, OpenBabel 3.1.1, ADMET-AI v2.1

**Environment Setup:**
```bash
conda env create -f environment.yml
conda activate timp2-docking
