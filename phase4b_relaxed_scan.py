import json
import shutil
import subprocess
from pathlib import Path
import numpy as np

HARTREE_TO_KCAL = 627.509474

W = Path("biohybrid_screen_v2")
W_MULTI = Path("relaxed_hop_v2_multistart")
W_MULTI.mkdir(exist_ok=True)
RES_DIR = Path("results")
RES_DIR.mkdir(exist_ok=True)

XTB_BIN = shutil.which("xtb")
if not XTB_BIN:
    raise RuntimeError("xtb executable not found in PATH")


def get_energy(out_text):
    for line in reversed(out_text.splitlines()):
        if "TOTAL ENERGY" in line:
            for tok in line.split():
                try:
                    val = float(tok)
                    if val < -1.0:
                        return val
                except ValueError:
                    continue
    raise ValueError("TOTAL ENERGY not found in xtb output")


def parse_xyz(xyz_path):
    lines = Path(xyz_path).read_text(encoding="utf-8").strip().splitlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2 : 2 + n]:
        sym, xs, ys, zs = line.split()[:4]
        atoms.append((sym, np.array([float(xs), float(ys), float(zs)])))
    return atoms


def write_xyz(path, atoms, comment=""):
    lines = [str(len(atoms)), comment]
    for sym, vec in atoms:
        x, y, z = vec
        lines.append(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


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
            shifted = shifted @ (np.eye(3) + np.sin(ang) * K + (1.0 - np.cos(ang)) * (K @ K)).T
    return [(sym, new_vec) for (sym, _), new_vec in zip(lig_atoms, shifted)]


def distal_anchor_indices(atoms_list, start_1based):
    o_coords = [vec for sym, vec in atoms_list if sym == "O"]
    anchors = []
    for idx, (sym, pos) in enumerate(atoms_list):
        if sym in ("C", "N"):
            min_do = min(float(np.linalg.norm(pos - oc)) for oc in o_coords)
            if min_do > 2.35:
                anchors.append(start_1based + idx)
    return anchors


def min_mg_o_distance(atoms):
    mg_pos = atoms[0][1]
    o_dists = [float(np.linalg.norm(pos - mg_pos)) for sym, pos in atoms[1:] if sym == "O"]
    return min(o_dists) if o_dists else 99.0


def build_host(use_bridge):
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
    return host, dme, anchors


def fresh_geometry(xm, host, dme, y_offset=-1.65, z_offset=0.0):
    dme_placed = [
        (s, x + np.array([xm, y_offset, z_offset]))
        for s, x in orient_ligand(dme, np.array([0.0, -1.0, 0.0]))
    ]
    return [("Mg", np.array([xm, 0.0, 0.0]))] + [(s, p.copy()) for s, p in host] + dme_placed


def shift_from_neighbor(xm, neighbor_atoms, n_host):
    prev_x = float(neighbor_atoms[0][1][0])
    dx = xm - prev_x
    shifted = []
    for idx, (s, pos) in enumerate(neighbor_atoms):
        if idx == 0:
            shifted.append(("Mg", np.array([xm, float(pos[1]), float(pos[2])])))
        elif idx >= 1 + n_host:
            shifted.append((s, pos + np.array([dx, 0.0, 0.0])))
        else:
            shifted.append((s, pos.copy()))
    return shifted


def run_xtb_opt(work_dir, atoms, anchors, xm, label, fix_mg=True):
    work_dir.mkdir(parents=True, exist_ok=True)
    write_xyz(work_dir / "input.xyz", atoms, f"{label} x_mg={xm:.3f}")
    fix_list = (["1"] if fix_mg else []) + [str(a) for a in anchors]
    fix_str = ",".join(fix_list)
    (work_dir / "xcontrol").write_text(f"$fix\n  atoms: {fix_str}\n$end\n", encoding="utf-8")
    res = subprocess.run(
        [XTB_BIN, "input.xyz", "--opt", "normal", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"],
        cwd=work_dir, capture_output=True, text=True
    )
    out = res.stdout + "\n" + res.stderr
    (work_dir / "xtb.out").write_text(out, encoding="utf-8")
    opt_xyz = work_dir / "xtbopt.xyz"
    if res.returncode != 0 or not opt_xyz.exists():
        return None, None
    opt_atoms = parse_xyz(opt_xyz)
    if min_mg_o_distance(opt_atoms) < 1.60:
        return None, None
    return get_energy(out), opt_atoms


def run_multistart_envelope(use_bridge=False):
    label = "hybrid_dma_envelope" if use_bridge else "pure_betaine_envelope"
    print(f"\n=== RUNNING MULTI-START ENVELOPE: {label} ===", flush=True)
    host, dme, anchors = build_host(use_bridge)
    n_host = len(host)
    xs = [round(float(v), 3) for v in np.linspace(-1.80, 1.80, 15)]

    fwd_geoms = {}
    prev_f = None
    for k, xm in enumerate(xs):
        init_g = fresh_geometry(xm, host, dme) if prev_f is None else shift_from_neighbor(xm, prev_f, n_host)
        e, opt_g = run_xtb_opt(W_MULTI / f"{label}_pt{k:02d}_s0_fwd", init_g, anchors, xm, "fwd")
        if opt_g is not None:
            prev_f = opt_g
            fwd_geoms[k] = (e, opt_g)

    rev_geoms = {}
    prev_r = None
    for k in reversed(range(len(xs))):
        xm = xs[k]
        init_g = fresh_geometry(xm, host, dme) if prev_r is None else shift_from_neighbor(xm, prev_r, n_host)
        e, opt_g = run_xtb_opt(W_MULTI / f"{label}_pt{k:02d}_s1_rev", init_g, anchors, xm, "rev")
        if opt_g is not None:
            prev_r = opt_g
            rev_geoms[k] = (e, opt_g)

    envelope_pts = []
    for k, xm in enumerate(xs):
        candidates = []
        if k in fwd_geoms:
            candidates.append(("fwd_chain", fwd_geoms[k][0], fwd_geoms[k][1]))
        if k in rev_geoms:
            candidates.append(("rev_chain", rev_geoms[k][0], rev_geoms[k][1]))

        seeds_pert = [(-1.75, 0.35, "pert_zpos"), (-1.75, -0.35, "pert_zneg")]
        for s_idx, (y_off, z_off, tag) in enumerate(seeds_pert, start=2):
            g_pert = fresh_geometry(xm, host, dme, y_offset=y_off, z_offset=z_off)
            e_p, opt_p = run_xtb_opt(W_MULTI / f"{label}_pt{k:02d}_s{s_idx}_{tag}", g_pert, anchors, xm, tag)
            if opt_p is not None:
                candidates.append((tag, e_p, opt_p))

        if not candidates:
            raise RuntimeError(f"All multi-start candidates failed at x_mg={xm}")

        best_src, best_e, best_geom = min(candidates, key=lambda item: item[1])
        spread_kcal = (max(c[1] for c in candidates) - best_e) * HARTREE_TO_KCAL
        write_xyz(W_MULTI / f"{label}_pt{k:02d}_BEST.xyz", best_geom, f"BEST ({best_src}) E={best_e:.8f}")

        envelope_pts.append({
            "x_A": xm,
            "abs_energy_eh": round(best_e, 8),
            "winning_seed": best_src,
            "seed_spread_kcal": round(spread_kcal, 2),
            "n_converged_seeds": len(candidates)
        })
        print(f"x={xm:+6.3f} A | E_min={best_e:.6f} Eh | seed={best_src:10s} | spread={spread_kcal:5.2f} kcal/mol", flush=True)

    e_min_global = min(p["abs_energy_eh"] for p in envelope_pts)
    for p in envelope_pts:
        p["rel_e_kcal"] = round((p["abs_energy_eh"] - e_min_global) * HARTREE_TO_KCAL, 3)

    left_end_dir = W_MULTI / f"{label}_opt_left_end"
    right_end_dir = W_MULTI / f"{label}_opt_right_end"
    left_atoms = parse_xyz(W_MULTI / f"{label}_pt00_BEST.xyz")
    right_atoms = parse_xyz(W_MULTI / f"{label}_pt{len(xs)-1:02d}_BEST.xyz")

    e_left_opt, g_left_opt = run_xtb_opt(left_end_dir, left_atoms, anchors, -1.80, "free_mg_left", fix_mg=False)
    e_right_opt, g_right_opt = run_xtb_opt(right_end_dir, right_atoms, anchors, 1.80, "free_mg_right", fix_mg=False)

    path_summary = {}
    if g_left_opt is not None and g_right_opt is not None:
        path_dir = W_MULTI / f"{label}_xtb_path"
        path_dir.mkdir(exist_ok=True)
        write_xyz(path_dir / "start.xyz", g_left_opt, "left endpoint")
        write_xyz(path_dir / "end.xyz", g_right_opt, "right endpoint")
        fix_host_only = ",".join(map(str, anchors))
        (path_dir / "xcontrol").write_text(f"$fix\n  atoms: {fix_host_only}\n$end\n", encoding="utf-8")
        res_p = subprocess.run(
            [XTB_BIN, "start.xyz", "--path", "end.xyz", "--chrg", "2", "--alpb", "thf", "--input", "xcontrol"],
            cwd=path_dir, capture_output=True, text=True
        )
        p_out = res_p.stdout + "\n" + res_p.stderr
        (path_dir / "xtb_path.out").write_text(p_out, encoding="utf-8")
        lx = round(float(g_left_opt[0][1][0]), 3)
        rx = round(float(g_right_opt[0][1][0]), 3)
        path_summary = {
            "left_min_eh": round(e_left_opt, 8),
            "left_min_mg_x_A": lx,
            "right_min_eh": round(e_right_opt, 8),
            "right_min_mg_x_A": rx,
            "xtb_path_returncode": res_p.returncode
        }
        print(f"Endpoint relaxation: Left Mg x={lx} A ({e_left_opt:.6f} Eh), Right Mg x={rx} A ({e_right_opt:.6f} Eh)", flush=True)

    return {
        "xs_valid": [p["x_A"] for p in envelope_pts],
        "abs_energies_eh": [p["abs_energy_eh"] for p in envelope_pts],
        "rel_e_kcal": [p["rel_e_kcal"] for p in envelope_pts],
        "points": envelope_pts,
        "endpoint_relaxation": path_summary
    }


def main():
    pure_env = run_multistart_envelope(use_bridge=False)
    hyb_env = run_multistart_envelope(use_bridge=True)

    out_file = RES_DIR / "multistart_hopping_envelope.json"
    out_file.write_text(json.dumps({
        "method": "GFN2-xTB / ALPB(THF), 4-seed multi-start lowest-energy envelope, Mg x-coordinate frozen ($fix)",
        "pure_betaine_envelope": pure_env,
        "hybrid_dma_envelope": hyb_env
    }, indent=2), encoding="utf-8")

    v2_path = RES_DIR / "verified_v2_data.json"
    if v2_path.exists():
        v2 = json.loads(v2_path.read_text(encoding="utf-8"))
        v2["multistart_hopping_envelope"] = {
            "pure_betaine_envelope": pure_env,
            "hybrid_dma_envelope": hyb_env
        }
        v2_path.write_text(json.dumps(v2, indent=2), encoding="utf-8")
    print(f"\nSaved multi-start envelope results to {out_file}")


if __name__ == "__main__":
    main()
