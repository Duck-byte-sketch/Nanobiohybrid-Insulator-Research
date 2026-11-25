import os
import shutil
import subprocess
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

HARTREE_TO_KCAL = 627.509474
WORK_DIR = Path("biohybrid_screen")
WORK_DIR.mkdir(exist_ok=True)

XTB_BIN = shutil.which("xtb")
assert XTB_BIN is not None, "xtb binary not found in PATH!"

# 1. Reference DME + four bio-hybrid insulator functional motifs
LIGANDS = {
    "DME":            {"smiles": "COCCOC",          "charge": 0,  "role": "Baseline Ether Solvent"},
    "Amide_NMA":      {"smiles": "CC(=O)NC",        "charge": 0,  "role": "Silk/Peptide C=O Donor"},
    "Alkoxy_OH":      {"smiles": "COCCO",           "charge": 0,  "role": "Cellulose/Chitosan O,O-Chelate"},
    "Carboxylate":    {"smiles": "CC(=O)[O-]",      "charge": -1, "role": "Alginate Anionic Anchor"},
    "Zwitterion_Bet": {"smiles": "C[N+](C)(C)CC(=O)[O-]", "charge": 0, "role": "Phospholipid/Betaine Dipole"},
}

def smiles_to_xyz(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    conf = mol.GetConformer()
    atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        atoms.append((atom.GetSymbol(), np.array([pos.x, pos.y, pos.z])))
    return atoms

def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment]
    for sym, xyz in atoms:
        lines.append(f"{sym:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}")
    Path(path).write_text("\n".join(lines) + "\n")

def run_xtb_opt(folder, atoms, charge, solvent="thf"):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    inp = folder / "input.xyz"
    write_xyz(inp, atoms, f"charge={charge}")
    cmd = [XTB_BIN, "input.xyz", "--opt", "normal", "--chrg", str(charge), "--alpb", solvent]
    res = subprocess.run(cmd, cwd=folder, capture_output=True, text=True)
    (folder / "xtb.out").write_text(res.stdout + "\n" + res.stderr)
    opt_xyz = folder / "xtbopt.xyz"
    if not opt_xyz.exists():
        raise RuntimeError(f"xTB failed in {folder}; check {folder}/xtb.out")
    # Parse energy from xtb.out
    energy = None
    for line in reversed(res.stdout.splitlines()):
        if "TOTAL ENERGY" in line:
            for tok in line.split():
                try:
                    val = float(tok)
                    if val < -1.0:
                        energy = val
                        break
                except ValueError:
                    continue
        if energy is not None:
            break
    # Parse optimized coordinates
    lines = opt_xyz.read_text().strip().splitlines()
    n = int(lines[0].strip())
    opt_atoms = []
    for line in lines[2:2+n]:
        p = line.split()
        opt_atoms.append((p[0], np.array([float(p[1]), float(p[2]), float(p[3])])))
    return energy, opt_atoms

def build_octahedral_cluster(dme_atoms, third_lig_atoms):
    """Places Mg2+ at origin, two DMEs in equatorial/axial chelate pockets, and 3rd ligand in open pocket."""
    cluster = [("Mg", np.array([0.0, 0.0, 0.0]))]

    def orient_ligand(lig_atoms, target_dir, offset_dist=2.35):
        o_indices = [i for i, (s, _) in enumerate(lig_atoms) if s == "O"]
        coords = np.array([xyz for _, xyz in lig_atoms])
        o_center = coords[o_indices].mean(axis=0)
        shifted = coords - o_center
        # Push oxygen center toward target_dir at offset_dist, keeping heavy atoms outward
        com = shifted.mean(axis=0)
        if np.linalg.norm(com) > 1e-3:
            v_out = com / np.linalg.norm(com)
            # Align v_out with target_dir
            u = target_dir / np.linalg.norm(target_dir)
            axis = np.cross(v_out, u)
            if np.linalg.norm(axis) > 1e-6:
                axis /= np.linalg.norm(axis)
                angle = np.arccos(np.clip(np.dot(v_out, u), -1.0, 1.0))
                K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
                R = np.eye(3) + np.sin(angle)*K + (1 - np.cos(angle))*(K @ K)
                shifted = shifted @ R.T
        placed = shifted + (target_dir / np.linalg.norm(target_dir)) * offset_dist
        return [(lig_atoms[i][0], placed[i]) for i in range(len(lig_atoms))]

    dirs = [
        np.array([1.0, 0.2, 0.0]),
        np.array([-0.5, 0.866, 0.2]),
        np.array([-0.5, -0.866, -0.2]),
    ]
    cluster.extend(orient_ligand(dme_atoms, dirs[0], 2.45))
    cluster.extend(orient_ligand(dme_atoms, dirs[1], 2.45))
    cluster.extend(orient_ligand(third_lig_atoms, dirs[2], 2.35))
    return cluster

def analyze_mg_o(atoms):
    mg_xyz = atoms[0][1]
    dists = sorted([np.linalg.norm(xyz - mg_xyz) for s, xyz in atoms[1:] if s == "O"])
    cn = sum(1 for d in dists if d < 2.45)
    return cn, [round(d, 2) for d in dists]

print("=== PHASE 2: NANOBIOHYBRID INSULATOR BINDING-SITE SCREEN (GFN2-xTB / ALPB-THF) ===", flush=True)

# Step 1: Optimize free ligands in ALPB(THF)
free_e = {}
free_geom = {}
for name, info in LIGANDS.items():
    atoms = smiles_to_xyz(info["smiles"])
    e, opt_atoms = run_xtb_opt(WORK_DIR / f"free_{name}", atoms, info["charge"])
    free_e[name] = e
    free_geom[name] = opt_atoms
    print(f"[Free Ligand] {name:15s} (q={info['charge']:+d}) : {e:14.6f} Eh", flush=True)

# Step 2: Optimize [Mg(DME)3]2+ baseline and [Mg(DME)2(BioMotif)] complexes
print("\n--- Optimizing [Mg(DME)2(Ligand)] Complexes in ALPB(THF) ---", flush=True)
complex_e = {}
results = []

for name, info in LIGANDS.items():
    init_cluster = build_octahedral_cluster(free_geom["DME"], free_geom[name])
    tot_charge = 2 + info["charge"]
    e_comp, opt_comp = run_xtb_opt(WORK_DIR / f"complex_Mg_DME2_{name}", init_cluster, tot_charge)
    complex_e[name] = e_comp
    cn, mg_o_dists = analyze_mg_o(opt_comp)
    results.append((name, info["role"], tot_charge, e_comp, cn, mg_o_dists))
    print(f"[Complex] Mg(DME)2({name:14s}) | q={tot_charge:+d} | CN={cn} | Mg-O (A)={mg_o_dists}", flush=True)

# Step 3: Compute DME-displacement energy: [Mg(DME)3]2+ + Motif -> [Mg(DME)2(Motif)] + DME
e_ref_complex = complex_e["DME"]
e_ref_dme = free_e["DME"]

print("\n=== SUMMARY: DME-DISPLACEMENT THERMODYNAMICS IN ALPB(THF) ===")
print(f"{'Motif':16s} | {'Bio-Analogue':28s} | {'CN':2s} | {'dE_swap (kcal/mol)':18s} | {'Mg-O Shell (A)'}")
print("-" * 105)
for name, role, q, e_comp, cn, dists in results:
    de_swap = ((e_comp + e_ref_dme) - (e_ref_complex + free_e[name])) * HARTREE_TO_KCAL
    print(f"{name:16s} | {role:28s} | {cn:2d} | {de_swap:18.2f} | {dists}")
