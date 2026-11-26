import json, shutil, subprocess
from pathlib import Path
import numpy as np

HARTREE_TO_KCAL = 627.509474
W = Path("biohybrid_screen_v2")
W_RED = Path("reduction_alpb_thf")
W_HOP = Path("relaxed_hop_v2")
W_HOP.mkdir(exist_ok=True)
RES_DIR = Path("results")
XTB_BIN = shutil.which("xtb")

data = json.loads((RES_DIR / "verified_v2_data.json").read_text())

# 1. Extract Mulliken spin populations (Mg vs. ligand) from doublet states
print("=== EXTRACTING DOUBLET SPIN LOCALIZATION (Mg vs. LIGAND) ===", flush=True)
for name, entry in data["vertical_reduction_alpb_thf"].items():
    safe = name.replace("[", "").replace("]", "").replace("+", "p")
    d = W_RED / f"{safe}_red"
    q0 = entry["charge_ox"]
    res = subprocess.run([XTB_BIN, "geom.xyz", "--sp", "--chrg", str(q0 - 1), "--uhf", "1", "--alpb", "thf", "--pop"],
                         cwd=d, capture_output=True, text=True)
    out = res.stdout + "\n" + res.stderr
    (d / "xtb_pop.out").write_text(out)
    spins = []
    in_table = False
    for line in out.splitlines():
        if "Mulliken/CM5 charges" in line:
            in_table = True
            continue
        if in_table:
            parts = line.split()
            if len(parts) >= 7 and parts.pop(0).isdigit():
                # Remaining parts: Z, sym, q_mull, q_cm5, covCN, alpha-beta_spin (or similar)
                z_str, sym, *nums = parts
                val = float(nums.pop())
                spins.append((sym, val))
            elif len(spins) > 0:
                break
    mg_sp = sum(v for s, v in spins if s == "Mg")
    lig_sp = sum(v for s, v in spins if s != "Mg")
    top = sorted(spins, key=lambda pair: abs(dict([pair]).popitem().pop()), reverse=True)[:4]
    sp_info = {
        "mg_spin": round(mg_sp, 3),
        "ligand_spin": round(lig_sp, 3),
        "top_atoms": [(s, round(v, 3)) for s, v in top],
    }
    entry["spin_localization"] = sp_info
    print(f"{name:27s} | Mg_spin={sp_info['mg_spin']:+6.3f} | Lig_spin={sp_info['ligand_spin']:+6.3f} | top={sp_info['top_atoms']}", flush=True)

# 2. Forward & Reverse Relaxed [Mg(DME)]2+ Hopping Scans + Clean Barrier Extraction
def get_energy(out_text):
    for line in reversed(out_text.splitlines()):
        if "TOTAL ENERGY" in line:
            for tok in line.split():
                try:
                    v = float(tok)
                    if v < -1.0:
                        return v
                except ValueError:
                    continue
    raise ValueError("TOTAL ENERGY not found")

def parse_xyz(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    _comment = lines.pop(0)
    atoms = []
    for l in lines[:n]:
        sym, xs, ys, zs = l.split()
        atoms.append((sym, np.array([float(xs), float(ys), float(zs)])))
    return atoms

def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment]
    for sym, vec in atoms:
        x, y, z = vec
        lines.append(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    Path(path).write_text("\n".join(lines) + "\n")

def orient_ligand(lig_atoms, target_dir):
    o_pts = np.array([vec for sym, vec in lig_atoms if sym == "O"])
    coords = np.array([vec for sym, vec in lig_atoms])
    shifted = coords - o_pts.mean(axis=0)
    com = shifted.mean(axis=0)
    u = target_dir / np.linalg.norm(target_dir)
    if np.linalg.norm(com) > 1e-3:
        v_out = com / np.linalg.norm(com)
        ax = np.cross(v_out, u)
        if np.linalg.norm(ax) > 1e-6:
            ax /= np.linalg.norm(ax)
            ax_x, ax_y, ax_z = ax
            ang = np.arccos(np.clip(np.dot(v_out, u), -1.0, 1.0))
            K = np.array([[0.0, -ax_z, ax_y], [ax_z, 0.0, -ax_x], [-ax_y, ax_x, 0.0]])
            shifted = shifted @ (np.eye(3) + np.sin(ang)*K + (1.0 - np.cos(ang))*(K @ K)).T
    out = []
    for (sym, _), new_vec in zip(lig_atoms, shifted):
        out.append((sym, new_vec))
    return out

def distal_anchor_indices(atoms_list, start_1based):
    o_coords = [vec for sym, vec in atoms_list if sym == "O"]
    anchors = []
    for idx, (sym, pos) in enumerate(atoms_list):
        if sym in ("C", "N"):
            min_do = min(float(np.linalg.norm(pos - oc)) for oc in o_coords)
            if min_do > 2.35:
                anchors.append(start_1based + idx)
    return anchors

def extract_barriers_clean(xs, energies_eh, min_r_mgo_list, clash_cutoff=1.65):
    vx_list, ve_list = [], []
    for x, e, r in zip(xs, energies_eh, min_r_mgo_list):
        if r >= clash_cutoff:
            vx_list.append(float(x))
            ve_list.append(float(e) * HARTREE_TO_KCAL)
    ve_arr = np.array(ve_list)
    ve_arr -= np.min(ve_arr)
    ve_clean = [float(v) for v in ve_arr]

    minima = []
    for i in range(1, len(ve_clean) - 1):
        prev_v = list(ve_clean[:i]).pop()
        curr_v = list(ve_clean[:i+1]).pop()
        next_v = list(ve_clean[:i+2]).pop()
        if curr_v <= prev_v and curr_v <= next_v:
            minima.append(i)
    if len(minima) < 2:
        mid = len(ve_clean) // 2
        l_idx = int(np.argmin(ve_arr[:mid]))
        r_idx = int(mid + np.argmin(ve_arr[mid:]))
    else:
        l_idx = int(list(minima[:1]).pop())
        r_idx = int(list(minima).pop())

    seg = ve_arr[l_idx:r_idx+1]
    ts_idx = int(l_idx + np.argmax(seg))
    ea_fwd = float(ve_arr[ts_idx] - ve_arr[l_idx])
    ea_rev = float(ve_arr[ts_idx] - ve_arr[r_idx])
    return {
        "xs_valid": [round(x, 3) for x in vx_list],
        "rel_e_kcal": [round(e, 3) for e in ve_clean],
        "min_L_A": round(float(np.array(vx_list)[l_idx]), 3),
        "ts_A": round(float(np.array(vx_list)[ts_idx]), 3),
        "min_R_A": round(float(np.array(vx_list)[r_idx]), 3),
        "Ea_fwd_kcal": round(ea_fwd, 2),
        "Ea_rev_kcal": round(ea_rev, 2),
        "Ea_fwd_eV": round(ea_fwd / 23.0605, 3),
        "Ea_rev_eV": round(ea_rev / 23.0605, 3),
    }

def run_directional_scan(use_bridge=False, reverse=False):
    direction = "rev" if reverse else "fwd"
    tag = f"{'hybrid_dma' if use_bridge else 'pure_bet'}_{direction}"
    bet = parse_xyz(W / "free_Zwitterion_Bet/xtbopt.xyz")
    dma = parse_xyz(W / "free_Amide_DMA/xtbopt.xyz")
    dme = parse_xyz(W / "free_DME/xtbopt.xyz")

    bet1 = [(s, x + np.array([-3.6, 0.0, 0.0])) for s, x in orient_ligand(bet, np.array([-1.0, 0.15, 0.0]))]
    bet2 = [(s, x + np.array([ 3.6, 0.0, 0.0])) for s, x in orient_ligand(bet, np.array([ 1.0, -0.15, 0.0]))]
    bridge = [(s, x + np.array([0.0, 2.35, 0.0])) for s, x in orient_ligand(dma, np.array([0.0, 1.0, 0.0]))] if use_bridge else []

    host = bet1 + bet2 + bridge
    anchors = distal_anchor_indices(bet1, 2) + distal_anchor_indices(bet2, 2 + len(bet1))
    if use_bridge:
        anchors += distal_anchor_indices(bridge, 2 + len(bet1) + len(bet2))

    xs = [float(v) for v in np.linspace(-1.80, 1.80, 15)]
    if reverse:
        xs.reverse()

    energies, min_mgo = [], []
    prev_coords = None
    for k, xm in enumerate(xs):
        d = W_HOP / f"{tag}_{k}"
        d.mkdir(exist_ok=True)
        if prev_coords is None:
            dme_placed = [(s, x + np.array([xm, -1.65, 0.0])) for s, x in orient_ligand(dme, np.array([0.0, -1.0, 0.0]))]
            sys_atoms = [("Mg", np.array([xm, 0.0, 0.0]))] + host + dme_placed
        else:
            sys_atoms = []
            first_sym, first_pos = list(prev_coords[:1]).pop()
            prev_x, _, _ = first_pos
            dx = xm - float(prev_x)
            for idx, (s, pos) in enumerate(prev_coords):
                if idx == 0:
                    _, py, pz = pos
                    sys_atoms.append(("Mg", np.array([xm, float(py), float(pz)])))
                elif idx >= 1 + len(host):
                    sys_atoms.append((s, pos + np.array([dx, 0.0, 0.0])))
                else:
                    sys_atoms.append((s, pos))

        write_xyz(d / "input.xyz", sys_atoms, f"x_mg={xm:.2f}")
        fix_str = "1," + ",".join(map(str, anchors))
        (d / "xcontrol").write_text(f"$fix\n  atoms: {fix_str}\n$end\n")
        res = subprocess.run([XTB_BIN, "input.xyz", "--opt", "normal", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"],
                             cwd=d, capture_output=True, text=True)
        out_txt = res.stdout + "\n" + res.stderr
        (d / "xtb.out").write_text(out_txt)
        e = get_energy(out_txt)
        opt_atoms = parse_xyz(d / "xtbopt.xyz")
        prev_coords = opt_atoms
        _, mg_p = list(opt_atoms[:1]).pop()
        r_min = min(float(np.linalg.norm(pos - mg_p)) for sym, pos in opt_atoms[1:] if sym == "O")
        energies.append(e)
        min_mgo.append(r_min)

    if reverse:
        xs.reverse()
        energies.reverse()
        min_mgo.reverse()

    return extract_barriers_clean(xs, energies, min_mgo)

print("\n=== PART 2: FORWARD & REVERSE RELAXED [Mg(DME)]2+ HOPPING SCANS ===", flush=True)
hop_results = {
    "pure_betaine_fwd": run_directional_scan(use_bridge=False, reverse=False),
    "pure_betaine_rev": run_directional_scan(use_bridge=False, reverse=True),
    "hybrid_dma_fwd":   run_directional_scan(use_bridge=True,  reverse=False),
    "hybrid_dma_rev":   run_directional_scan(use_bridge=True,  reverse=True),
}

for k, v in hop_results.items():
    print(f"{k:20s} | Min_L={v['min_L_A']:+.2f} A, TS={v['ts_A']:+.2f} A, Min_R={v['min_R_A']:+.2f} A | "
          f"Ea(fwd)={v['Ea_fwd_kcal']:5.2f} kcal/mol ({v['Ea_fwd_eV']:.2f} eV) | "
          f"Ea(rev)={v['Ea_rev_kcal']:5.2f} kcal/mol ({v['Ea_rev_eV']:.2f} eV)", flush=True)

data["relaxed_hopping_scans"] = hop_results
(RES_DIR / "verified_v2_data.json").write_text(json.dumps(data, indent=2) + "\n")
print("\nAll verified results saved to results/verified_v2_data.json!", flush=True)
