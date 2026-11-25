import os, shutil, subprocess
from pathlib import Path
import numpy as np
from pyscf import gto, dft

os.environ["PYSCF_TMPDIR"] = os.path.expanduser("~/pyscf_scratch")
HARTREE_TO_KCAL = 627.509474
HARTREE_TO_EV = 27.211386245988
W = Path("biohybrid_screen_v2")
W_HOP = Path("relaxed_hop_v2")
W_HOP.mkdir(exist_ok=True)
XTB_BIN = shutil.which("xtb")

def get_xtb_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            for tok in line.split():
                try:
                    v = float(tok)
                    if v < -1.0:
                        return v
                except ValueError:
                    continue
    raise ValueError(f"Could not parse TOTAL ENERGY in {out_path}")

def parse_xyz(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.strip())
    return [(l.split(), np.array([float(l.split()), float(l.split()), float(l.split())])) for l in lines[2:2+n]]

def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment] + [f"{s:2s} {x:12.6f} {x:12.6f} {x:12.6f}" for s, x in atoms]
    Path(path).write_text("\n".join(lines) + "\n")

def orient_ligand(lig_atoms, target_dir, offset_dist=2.35):
    o_idx = [i for i, (s, _) in enumerate(lig_atoms) if s == "O"]
    coords = np.array([x for _, x in lig_atoms])
    shifted = coords - coords[o_idx].mean(axis=0)
    com = shifted.mean(axis=0)
    u = target_dir / np.linalg.norm(target_dir)
    if np.linalg.norm(com) > 1e-3:
        v_out = com / np.linalg.norm(com)
        ax = np.cross(v_out, u)
        if np.linalg.norm(ax) > 1e-6:
            ax /= np.linalg.norm(ax)
            ang = np.arccos(np.clip(np.dot(v_out, u), -1.0, 1.0))
            K = np.array([[0, -ax, ax], [ax, 0, -ax], [-ax, ax, 0]])
            shifted = shifted @ (np.eye(3) + np.sin(ang)*K + (1 - np.cos(ang))*(K @ K)).T
    return [(lig_atoms[i], shifted[i] + u * offset_dist) for i in range(len(lig_atoms))]

def run_relaxed_hop(use_bridge=False):
    tag = "relaxed_hybrid_dma" if use_bridge else "relaxed_pure_bet"
    bet = parse_xyz(W / "free_Zwitterion_Bet/xtbopt.xyz")
    dma = parse_xyz(W / "free_Amide_DMA/xtbopt.xyz")
    dme = parse_xyz(W / "free_DME/xtbopt.xyz")

    bet1 = [(s, x + np.array([-3.4, 0.0, 0.0])) for s, x in orient_ligand(bet, np.array([-1.0, 0.2, 0.0]), 0.0)]
    bet2 = [(s, x + np.array([ 3.4, 0.0, 0.0])) for s, x in orient_ligand(bet, np.array([ 1.0, -0.2, 0.0]), 0.0)]
    bridge = [(s, x + np.array([0.0, 2.2, 0.0])) for s, x in orient_ligand(dma, np.array([0.0, 1.0, 0.0]), 0.0)] if use_bridge else []

    host = bet1 + bet2 + bridge
    anchor_indices = [i + 2 for i, (s, _) in enumerate(host) if s in ("C", "N")]

    xs = np.linspace(-1.45, 1.45, 13)
    energies = []
    for k, xm in enumerate(xs):
        d = W_HOP / f"{tag}_{k}"
        d.mkdir(exist_ok=True)
        dme_placed = [(s, x + np.array([xm, -1.55, 0.0])) for s, x in orient_ligand(dme, np.array([0.0, -1.0, 0.0]), 0.0)]
        sys_atoms = [("Mg", np.array([xm, 0.0, 0.0]))] + host + dme_placed
        write_xyz(d / "input.xyz", sys_atoms, f"x_mg={xm:.2f}")
        fix_str = "1," + ",".join(map(str, anchor_indices))
        (d / "xcontrol").write_text(f"$fix\n  atoms: {fix_str}\n$end\n")
        res = subprocess.run([XTB_BIN, "input.xyz", "--opt", "normal", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"],
                             cwd=d, capture_output=True, text=True)
        (d / "xtb.out").write_text(res.stdout + "\n" + res.stderr)
        energies.append(get_xtb_energy(d / "xtb.out"))

    e_kcal = np.array(energies) * HARTREE_TO_KCAL
    mid = len(xs) // 2
    l_idx = int(np.argmin(e_kcal[:mid]))
    r_idx = int(mid + 1 + np.argmin(e_kcal[mid+1:]))
    ts_idx = int(l_idx + np.argmax(e_kcal[l_idx:r_idx+1]))
    e_rel = e_kcal - e_kcal[l_idx]
    ea_fwd = float(e_rel[ts_idx] - e_rel[l_idx])
    ea_rev = float(e_rel[ts_idx] - e_rel[r_idx])
    return xs, e_rel, (float(xs[l_idx]), float(xs[ts_idx]), float(xs[r_idx]), ea_fwd, ea_rev)

print("=== STEP 3: RELAXED [Mg(DME)]2+ MIGRATION PROFILE ===", flush=True)
xs, ep, sp = run_relaxed_hop(False)
_,  eh, sh = run_relaxed_hop(True)
h0, h1, h2 = "x_Mg (A)", "Pure Betaine (kcal/mol)", "Betaine+DMA Bridge (kcal/mol)"
print(f"{h0:>8s} | {h1:>24s} | {h2:>29s}", flush=True)
print("-" * 68, flush=True)
for x, a, b in zip(xs, ep, eh):
    print(f"{x:8.2f} | {a:24.2f} | {b:29.2f}", flush=True)
print(f"\nPure Betaine   | Min_L={sp:+.2f} A, TS={sp:+.2f} A, Min_R={sp:+.2f} A | E_a(fwd)={sp:.2f} kcal/mol ({sp/23.06:.2f} eV)", flush=True)
print(f"Betaine + DMA  | Min_L={sh:+.2f} A, TS={sh:+.2f} A, Min_R={sh:+.2f} A | E_a(fwd)={sh:.2f} kcal/mol ({sh/23.06:.2f} eV)", flush=True)
print(f"Bridge Effect  | dE_a = {sh - sp:+.2f} kcal/mol ({(sh - sp)/23.06:+.2f} eV)", flush=True)

def xyz_to_atoms_str(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.strip())
    return "\n".join(lines[2:2+n])

def run_svp_pcm(xyz_path, charge, spin, label):
    mol = gto.M(atom=xyz_to_atoms_str(xyz_path), basis="def2-svp", charge=charge, spin=spin, verbose=0)
    mf = (dft.RKS(mol) if spin == 0 else dft.UKS(mol)).density_fit(auxbasis="def2-universal-jkfit")
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf.max_memory = 11000
    mf.conv_tol = 1e-6
    mf = mf.PCM()
    mf.with_solvent.method = "IEF-PCM"
    mf.with_solvent.eps = 7.58
    e = mf.kernel()
    print(f"[DFT-SVP-PCM(THF)] {label:25s} (q={charge:+d}, S={spin}) : {e:16.8f} Eh | conv={mf.converged}", flush=True)
    return e

print("\n=== FAST VERTICAL ELECTRON REDUCTION IN PCM(THF) (B3LYP-D3BJ/def2-SVP/IEF-PCM) ===", flush=True)
for tag, p in [("[Mg(DME)3]", "complex_Mg_DME2_DME/xtbopt.xyz"), ("[Mg(DME)2(DMA)]", "complex_Mg_DME2_Amide_DMA/xtbopt.xyz"), ("[Mg(DME)2(Bet)]", "complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz")]:
    e2 = run_svp_pcm(W / p, 2, 0, f"{tag}2+")
    e1 = run_svp_pcm(W / p, 1, 1, f"{tag}.+")
    vea = (e2 - e1) * HARTREE_TO_EV
    print(f"  -> VEA {tag:18s} = {vea:6.2f} eV", flush=True)
