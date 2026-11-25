import shutil, subprocess
from pathlib import Path
import numpy as np

HARTREE_TO_KCAL = 627.509474
W_PREV = Path("biohybrid_screen")
W4 = Path("phase4_transport")
XTB_BIN = shutil.which("xtb")

def parse_xyz(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    return [(l.split()[0], np.array([float(l.split()[1]), float(l.split()[2]), float(l.split()[3])])) for l in lines[:n]]

def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment] + [f"{s:2s} {x[0]:12.6f} {x[1]:12.6f} {x[2]:12.6f}" for s, x in atoms]
    Path(path).write_text("\n".join(lines) + "\n")

def align_site(lig_atoms, target_o, out_vec):
    o_idx = [i for i, (s, _) in enumerate(lig_atoms) if s == "O"]
    coords = np.array([x for _, x in lig_atoms])
    shifted = coords - coords[o_idx].mean(axis=0)
    com = shifted.mean(axis=0)
    if np.linalg.norm(com) > 1e-4:
        v1, v2 = com / np.linalg.norm(com), out_vec / np.linalg.norm(out_vec)
        ax = np.cross(v1, v2)
        if np.linalg.norm(ax) > 1e-6:
            ax /= np.linalg.norm(ax)
            ang = np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
            K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
            shifted = shifted @ (np.eye(3) + np.sin(ang)*K + (1 - np.cos(ang))*(K @ K)).T
    return [(lig_atoms[i][0], shifted[i] + target_o) for i in range(len(lig_atoms))]

def scan_wide(use_bridge=False):
    tag = "wide_hybrid" if use_bridge else "wide_pure"
    bet = parse_xyz(W_PREV / "free_Zwitterion_Bet/xtbopt.xyz")
    nma = parse_xyz(W_PREV / "free_Amide_NMA/xtbopt.xyz")
    # 7.4 A spacing between Betaine O-centers (x = -3.7 A and +3.7 A -> equilibrium Mg-O ~ 1.95 A at x = +/-1.75 A)
    host = align_site(bet, np.array([-3.7, 0.0, 0.0]), np.array([-1.0, 0.2, 0.0])) + \
           align_site(bet, np.array([ 3.7, 0.0, 0.0]), np.array([ 1.0, -0.2, 0.0]))
    if use_bridge:
        host += align_site(nma, np.array([0.0, 1.95, 0.0]), np.array([0.0, 1.0, 0.0]))
    heavy = [i + 2 for i, (s, _) in enumerate(host) if s != "H"]
    xs = np.linspace(-1.75, 1.75, 8)
    es = []
    for k, xm in enumerate(xs):
        d = W4 / f"{tag}_{k}"
        d.mkdir(exist_ok=True)
        write_xyz(d / "input.xyz", [("Mg", np.array([xm, 0.0, 0.0]))] + host)
        (d / "xcontrol").write_text(f"$fix\n  atoms: 1,{','.join(map(str, heavy))}\n$end\n")
        res = subprocess.run([XTB_BIN, "input.xyz", "--opt", "loose", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"],
                             cwd=d, capture_output=True, text=True)
        for line in reversed(res.stdout.splitlines()):
            if "TOTAL ENERGY" in line:
                es.append(float([t for t in line.split() if t.replace('.','',1).replace('-','',1).isdigit() and float(t) < -1.0][0]))
                break
    rel = (np.array(es) - min(es)) * HARTREE_TO_KCAL
    return xs, rel

xs, ep = scan_wide(False)
_,  eh = scan_wide(True)
print("=== WIDE POLYMER SPACING (7.4 A): Mg2+ HOPPING BARRIER (GFN2-xTB / ALPB-THF) ===")
print(f"{'x_Mg (A)':>8s} | {'Pure Betaine (kcal/mol)':>24s} | {'Betaine+Amide Bridge (kcal/mol)':>31s}")
print("-" * 70)
for x, a, b in zip(xs, ep, eh):
    print(f"{x:8.2f} | {a:24.2f} | {b:31.2f}")
print(f"\nWide-Spacing Barrier (Pure Betaine)         : {np.max(ep):6.2f} kcal/mol")
print(f"Wide-Spacing Barrier (Betaine + Amide)      : {np.max(eh):6.2f} kcal/mol")
print(f"Cooperative Amide Bridge Barrier Reduction  : {np.max(eh) - np.max(ep):+6.2f} kcal/mol")
