import math, subprocess, shutil, sys
from pathlib import Path

HARTREE_TO_KCAL = 627.509474
BASE = Path("mg_benchmark")

xtb_bin = shutil.which("xtb")
if not xtb_bin:
    sys.exit("ERROR: 'xtb' not found on PATH. Run 'conda activate mg_bench' first!")
print(f"Using xtb binary: {xtb_bin}")

def write_xyz(path, elems, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

# Ideal gauche bidentate DME coordinates (16 atoms) oriented to chelate Mg at (0,0,0) along +X
# Both ether oxygens sit at x = +1.65, y = +/- 1.30 (Mg-O = 2.10 A, O-Mg-O bite angle = 76 deg)
gauche_dme_elems = [
    "O", "O", "C", "C", "C", "C",
    "H", "H", "H", "H", "H", "H", "H", "H", "H", "H"
]
gauche_dme_coords = [
    (1.65,  1.30,  0.00),  # O1
    (1.65, -1.30,  0.00),  # O2
    (2.85,  0.72,  0.35),  # C2 (ethylene bridge)
    (2.85, -0.72, -0.35),  # C3 (ethylene bridge)
    (1.55,  2.72, -0.15),  # C1 (methyl)
    (1.55, -2.72,  0.15),  # C4 (methyl)
    (3.72,  1.25,  0.00),  # H on C2
    (2.90,  0.75,  1.45),  # H on C2
    (3.72, -1.25,  0.00),  # H on C3
    (2.90, -0.75, -1.45),  # H on C3
    (0.62,  3.05,  0.32),  # H on C1
    (1.55,  3.00, -1.21),  # H on C1
    (2.41,  3.20,  0.34),  # H on C1
    (0.62, -3.05, -0.32),  # H on C4
    (1.55, -3.00,  1.21),  # H on C4
    (2.41, -3.20, -0.34),  # H on C4
]

def rot_euler(coords, rx, ry, rz):
    ax, ay, az = map(math.radians, (rx, ry, rz))
    cx, sx = math.cos(ax), math.sin(ax)
    cy, sy = math.cos(ay), math.sin(ay)
    cz, sz = math.cos(az), math.sin(az)
    out = []
    for x, y, z in coords:
        # Rx
        y1, z1 = y*cx - z*sx, y*sx + z*cx
        # Ry
        x2, z2 = x*cy + z1*sy, -x*sy + z1*cy
        # Rz
        x3, y3 = x2*cz - y1*sz, x2*sz + y1*cz
        out.append((x3, y3, z2))
    return out

# 3 orthogonal propeller-like chelate rings (saturates octahedral CN=6) + 3 second-shell DMEs
dme_slots = [
    ((0,   0,   0), (0.0,  0.0,  0.0)),   # Shell 1: XY plane
    ((90,  0, 120), (0.0,  0.0,  0.0)),   # Shell 1: tilted 120 deg
    ((-90, 0, 240), (0.0,  0.0,  0.0)),   # Shell 1: tilted 240 deg
    ((0,  90,   0), (0.0,  0.0,  4.8)),   # Shell 2: +Z
    ((0, -90,   0), (0.0,  0.0, -4.8)),   # Shell 2: -Z
    ((0,   0,  60), (2.4,  4.2,  0.0)),   # Shell 2: outer equatorial
]

for n in range(4, 7):
    elems, coords = ["Mg"], [(0.0, 0.0, 0.0)]
    for i in range(n):
        angs, (sx, sy, sz) = dme_slots[i]
        rc = rot_euler(gauche_dme_coords, *angs)
        elems.extend(gauche_dme_elems)
        coords.extend([(x + sx, y + sy, z + sz) for x, y, z in rc])
    write_xyz(BASE / f"mg_dme{n}_chelate.xyz", elems, coords, f"Octahedral chelated [Mg(DME)_{n}]2+")

def run_xtb(work_dir, in_file, chrg, solv_flag):
    shutil.copy(in_file, work_dir / "input.xyz")
    cmd = [xtb_bin, "input.xyz", "--opt", "tight", "--gfn", "2", "--chrg", str(chrg)] + solv_flag
    res = subprocess.run(cmd, cwd=work_dir, capture_output=True, text=True)
    log = res.stdout + "\n" + res.stderr
    (work_dir / "final.out").write_text(log)
    if res.returncode != 0 or not (work_dir / "xtbopt.xyz").exists():
        print(f"[ERROR] xtb failed in {work_dir}:\n" + "\n".join(log.splitlines()[-15:]))
        sys.exit(1)
    shutil.copy(work_dir / "xtbopt.xyz", work_dir / "final.xyz")

# Re-optimize DME clusters (n=4,5,6) and Urea clusters (n=4,5,6) cleanly
for n in range(4, 7):
    # DME from new gauche chelate guess
    run_xtb(BASE / f"mg_dme{n}_solv", BASE / f"mg_dme{n}_chelate.xyz", 2, ["--alpb", "thf"])
    run_xtb(BASE / f"mg_dme{n}_gas",  BASE / f"mg_dme{n}_solv/final.xyz", 2, [])
    # Urea from already-verified fixed guess
    run_xtb(BASE / f"mg_urea{n}_solv", BASE / f"mg_urea{n}_fixed.xyz", 2, ["--alpb", "dmso"])
    run_xtb(BASE / f"mg_urea{n}_gas",  BASE / f"mg_urea{n}_solv/final.xyz", 2, [])

def read_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines.strip())
    elems, coords = [], []
    for line in lines[2:2+n]:
        p = line.split()
        elems.append(p)
        coords.append((float(p), float(p), float(p)))
    return elems, coords

def get_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            return float(line.split())
    raise ValueError(f"No TOTAL ENERGY in {out_path}")

print("\n=== 1. Verified Coordination Numbers (Mg-O < 2.4 A) ===")
for lig in ("dme", "urea"):
    for n in range(4, 7):
        for phase in ("solv", "gas"):
            elems, coords = read_xyz(BASE / f"mg_{lig}{n}_{phase}" / "final.xyz")
            mx, my, mz = coords
            o_dists = [
                math.sqrt((cx-mx)**2 + (cy-my)**2 + (cz-mz)**2)
                for e, (cx, cy, cz) in zip(elems, coords) if e == "O"
            ]
            cn = sum(d < 2.40 for d in o_dists)
            near = sorted(round(d, 2) for d in o_dists)[:8]
            print(f"mg_{lig}{n}_{phase:4s} | CN(Mg-O)={cn} | nearest O (A): {near}")

print("\n=== 2. GFN2-xTB Stepwise Desolvation Energies (kcal/mol) ===")
for phase in ("solv", "gas"):
    print(f"\n--- Phase: {phase.upper()} ---")
    for lig in ("dme", "urea"):
        e_lig = get_energy(BASE / f"{lig}_{phase}" / "final.out")
        e4 = get_energy(BASE / f"mg_{lig}4_{phase}" / "final.out")
        e5 = get_energy(BASE / f"mg_{lig}5_{phase}" / "final.out")
        e6 = get_energy(BASE / f"mg_{lig}6_{phase}" / "final.out")
        de_6_5 = ((e5 + e_lig) - e6) * HARTREE_TO_KCAL
        de_5_4 = ((e4 + e_lig) - e5) * HARTREE_TO_KCAL
        print(f"{lig.upper():4s} | n=6->5: {de_6_5:8.2f} kcal/mol | n=5->4: {de_5_4:8.2f} kcal/mol")
