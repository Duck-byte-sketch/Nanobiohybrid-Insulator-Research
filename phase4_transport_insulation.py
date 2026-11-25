import shutil
import subprocess
from pathlib import Path
import numpy as np
from pyscf import gto, dft

HARTREE_TO_EV = 27.211386245988
HARTREE_TO_KCAL = 627.509474
W_PREV = Path("biohybrid_screen")
W4 = Path("phase4_transport")
W4.mkdir(exist_ok=True)

XTB_BIN = shutil.which("xtb")
assert XTB_BIN is not None, "xtb not found"

def xyz_to_atoms(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    return "\n".join(lines[:n])

def parse_xyz(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    atoms = []
    for line in lines[:n]:
        p = line.split()
        atoms.append((p[0], np.array([float(p[1]), float(p[2]), float(p[3])])))
    return atoms

def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment]
    for s, xyz in atoms:
        lines.append(f"{s:2s} {xyz[0]:12.6f} {xyz[1]:12.6f} {xyz[2]:12.6f}")
    Path(path).write_text("\n".join(lines) + "\n")

def compute_homo_lumo(xyz_path, charge, label):
    mol = gto.M(atom=xyz_to_atoms(xyz_path), basis="def2-svp", charge=charge, spin=0, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf.max_memory = 700
    mf.direct_scf = True
    mf.direct_scf_tol = 1e-10
    mf.grids.level = 2
    mf.conv_tol = 1e-6
    mf.kernel()
    occ = mf.mo_occ > 0
    homo = mf.mo_energy[occ][-1] * HARTREE_TO_EV
    lumo = mf.mo_energy[~occ][0] * HARTREE_TO_EV
    gap = lumo - homo
    print(f"[B3LYP-D3BJ/def2-SVP] {label:22s} | HOMO={homo:7.2f} eV | LUMO={lumo:7.2f} eV | Gap={gap:6.2f} eV", flush=True)
    return homo, lumo, gap

print("=== PART A: ELECTRONIC INSULATION (B3LYP-D3BJ/def2-SVP HOMO-LUMO GAPS) ===", flush=True)
compute_homo_lumo(W_PREV / "free_DME/xtbopt.xyz", 0, "Free DME (Electrolyte)")
compute_homo_lumo(W_PREV / "free_Amide_NMA/xtbopt.xyz", 0, "Free Amide (NMA)")
compute_homo_lumo(W_PREV / "free_Zwitterion_Bet/xtbopt.xyz", 0, "Free Zwitterion_Bet")
compute_homo_lumo(W_PREV / "complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2, "[Mg(DME)2(Bet)]2+")

def align_ligand_to_site(lig_atoms, target_o_pos, outward_vec):
    o_idx = [i for i, (s, _) in enumerate(lig_atoms) if s == "O"]
    coords = np.array([xyz for _, xyz in lig_atoms])
    o_center = coords[o_idx].mean(axis=0)
    shifted = coords - o_center
    com = shifted.mean(axis=0)
    if np.linalg.norm(com) > 1e-4:
        v1 = com / np.linalg.norm(com)
        v2 = outward_vec / np.linalg.norm(outward_vec)
        axis = np.cross(v1, v2)
        if np.linalg.norm(axis) > 1e-6:
            axis /= np.linalg.norm(axis)
            ang = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
            K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            R = np.eye(3) + np.sin(ang)*K + (1 - np.cos(ang))*(K @ K)
            shifted = shifted @ R.T
    placed = shifted + target_o_pos
    return [(lig_atoms[i][0], placed[i]) for i in range(len(lig_atoms))]

def run_hopping_scan(use_amide_bridge=False):
    tag = "hybrid_betaine_amide" if use_amide_bridge else "pure_betaine"
    bet_atoms = parse_xyz(W_PREV / "free_Zwitterion_Bet/xtbopt.xyz")
    nma_atoms = parse_xyz(W_PREV / "free_Amide_NMA/xtbopt.xyz")

    # Place Betaine 1 oxygen center at x = -2.6 A, Betaine 2 oxygen center at x = +2.6 A
    bet1 = align_ligand_to_site(bet_atoms, np.array([-2.6, 0.0, 0.0]), np.array([-1.0, 0.3, 0.0]))
    bet2 = align_ligand_to_site(bet_atoms, np.array([ 2.6, 0.0, 0.0]), np.array([ 1.0, -0.3, 0.0]))
    bridge = align_ligand_to_site(nma_atoms, np.array([0.0, 1.85, 0.0]), np.array([0.0, 1.0, 0.0])) if use_amide_bridge else []

    # Fix all heavy atoms of the polymer pocket anchors so only Mg2+ + local H/sidechains relax along x
    host_atoms = bet1 + bet2 + bridge
    heavy_indices = [i + 2 for i, (s, _) in enumerate(host_atoms) if s != "H"] # 1-indexed after Mg at index 1

    xs = np.linspace(-1.7, 1.7, 9)
    energies = []
    for step, x_mg in enumerate(xs):
        step_dir = W4 / f"{tag}_step{step}"
        step_dir.mkdir(parents=True, exist_ok=True)
        sys_atoms = [("Mg", np.array([x_mg, 0.0, 0.0]))] + host_atoms
        write_xyz(step_dir / "input.xyz", sys_atoms, f"x_mg={x_mg:.2f}")

        # Constrain Mg2+ (atom 1) and anchor heavy atoms so Mg2+ samples the exact migration coordinate
        fix_list = "1," + ",".join(map(str, heavy_indices))
        xcontrol = f"$fix\n  atoms: {fix_list}\n$end\n"
        (step_dir / "xcontrol").write_text(xcontrol)

        cmd = [XTB_BIN, "input.xyz", "--opt", "loose", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"]
        res = subprocess.run(cmd, cwd=step_dir, capture_output=True, text=True)
        e_val = None
        for line in reversed(res.stdout.splitlines()):
            if "TOTAL ENERGY" in line:
                for tok in line.split():
                    try:
                        v = float(tok)
                        if v < -1.0:
                            e_val = v
                            break
                    except ValueError:
                        continue
            if e_val is not None:
                break
        energies.append(e_val)

    e_kcal = (np.array(energies) - min(energies)) * HARTREE_TO_KCAL
    barrier = float(np.max(e_kcal))
    return xs, e_kcal, barrier

print("\n=== PART B: Mg2+ MIGRATION BARRIER SCAN (GFN2-xTB / ALPB-THF) ===", flush=True)
xs_pure, e_pure, bar_pure = run_hopping_scan(use_amide_bridge=False)
xs_hyb,  e_hyb,  bar_hyb  = run_hopping_scan(use_amide_bridge=True)

print(f"{'x_Mg (A)':>8s} | {'Pure Betaine (kcal/mol)':>24s} | {'Betaine+Amide Bridge (kcal/mol)':>31s}")
print("-" * 70)
for x, ep, eh in zip(xs_pure, e_pure, e_hyb):
    print(f"{x:8.2f} | {ep:24.2f} | {eh:31.2f}")

print(f"\n[BARRIER SUMMARY] Pure Zwitterionic Channel E_a   : {bar_pure:6.2f} kcal/mol")
print(f"[BARRIER SUMMARY] Hybrid Betaine-Amide Bridge E_a : {bar_hyb:6.2f} kcal/mol")
print(f"[BARRIER REDUCTION] Cooperative Amide Bridge dE_a : {bar_hyb - bar_pure:+6.2f} kcal/mol")
