import math, subprocess, shutil
from pathlib import Path

HARTREE_TO_KCAL = 627.509474
BASE = Path("mg_benchmark")

def read_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0].strip())
    elems, coords = [], []
    for line in lines[2:2+n]:
        p = line.split()
        elems.append(p[0])
        coords.append([float(p[1]), float(p[2]), float(p[3])])
    return elems, coords

def write_xyz(path, elems, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

def orient_dme(elems, coords):
    o_idx = [i for i, e in enumerate(elems) if e == "O"]
    c_idx = [i for i, e in enumerate(elems) if e != "O"]
    ox = sum(coords[i][0] for i in o_idx) / 2.0
    oy = sum(coords[i][1] for i in o_idx) / 2.0
    oz = sum(coords[i][2] for i in o_idx) / 2.0
    shifted = [[x - ox, y - oy, z - oz] for x, y, z in coords]
    cx = sum(shifted[i][0] for i in c_idx) / len(c_idx)
    if cx < 0:
        shifted = [[-x, y, z] for x, y, z in shifted]
    return shifted

def orient_urea(elems, coords):
    o_i = elems.index("O")
    c_i = elems.index("C")
    ox, oy, oz = coords[o_i]
    shifted = [[x - ox, y - oy, z - oz] for x, y, z in coords]
    vx, vy, vz = shifted[c_i]
    norm = math.sqrt(vx*vx + vy*vy + vz*vz)
    vx, vy, vz = vx/norm, vy/norm, vz/norm
    return [[x*vx + y*vy + z*vz, -x*vy + y*vx, z] for x, y, z in shifted]

def rot_xyz(coords, axis, angle_deg):
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)
    if axis == "z":
        return [[x*ca - y*sa, x*sa + y*ca, z] for x, y, z in coords]
    elif axis == "y":
        return [[x*ca + z*sa, y, -x*sa + z*ca] for x, y, z in coords]
    return [[x, y*ca - z*sa, y*sa + z*ca] for x, y, z in coords]

dme_e, dme_c = read_xyz(BASE / "dme_solv" / "final.xyz")
urea_e, urea_c = read_xyz(BASE / "urea_solv" / "final.xyz")
dme_o = orient_dme(dme_e, dme_c)
urea_o = orient_urea(urea_e, urea_c)

chelate_poses = [
    ("z", 0,   ( 2.15,  0.00,  0.00)),
    ("z", 120, (-1.08,  1.86,  0.00)),
    ("z", 240, (-1.08, -1.86,  0.00)),
    ("y", 90,  ( 0.00,  0.00,  5.80)),
    ("y", -90, ( 0.00,  0.00, -5.80)),
    ("z", 60,  ( 2.90,  5.00,  0.00)),
]

for n in range(4, 7):
    elems, coords = ["Mg"], [[0.0, 0.0, 0.0]]
    for i in range(n):
        ax, ang, (sx, sy, sz) = chelate_poses[i]
        rc = rot_xyz(dme_o, ax, ang)
        elems.extend(dme_e)
        coords.extend([[x + sx, y + sy, z + sz] for x, y, z in rc])
    write_xyz(BASE / f"mg_dme{n}_fixed.xyz", elems, coords, f"Chelated [Mg(DME)_{n}]2+")

urea_poses = [
    ("z", 0,   ( 2.05,  0.00,  0.00)),
    ("z", 180, (-2.05,  0.00,  0.00)),
    ("z", 90,  ( 0.00,  2.05,  0.00)),
    ("z", -90, ( 0.00, -2.05,  0.00)),
    ("y", -90, ( 0.00,  0.00,  2.05)),
    ("y", 90,  ( 0.00,  0.00, -2.05)),
]

for n in range(4, 7):
    elems, coords = ["Mg"], [[0.0, 0.0, 0.0]]
    for i in range(n):
        ax, ang, (sx, sy, sz) = urea_poses[i]
        rc = rot_xyz(urea_o, ax, ang)
        elems.extend(urea_e)
        coords.extend([[x + sx, y + sy, z + sz] for x, y, z in rc])
    write_xyz(BASE / f"mg_urea{n}_fixed.xyz", elems, coords, f"Coordinated [Mg(urea)_{n}]2+")

def opt_pair(lig, n, solv):
    in_xyz = BASE / f"mg_{lig}{n}_fixed.xyz"
    for phase, flag in (("solv", ["--alpb", solv]), ("gas", [])):
        d = BASE / f"mg_{lig}{n}_{phase}"
        shutil.copy(in_xyz, d / "input.xyz")
        cmd = ["xtb", "input.xyz", "--opt", "tight", "--gfn", "2", "--chrg", "2"] + flag
        res = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
        (d / "final.out").write_text(res.stdout + res.stderr)
        shutil.copy(d / "xtbopt.xyz", d / "final.xyz")
        in_xyz = d / "final.xyz"

for n in range(4, 7):
    opt_pair("dme", n, "thf")
    opt_pair("urea", n, "dmso")

def get_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            return float(line.split()[3])
    raise ValueError(f"No TOTAL ENERGY in {out_path}")

print("\n=== 1. Re-Audited Coordination Numbers (Mg-O < 2.4 A) ===")
for lig in ("dme", "urea"):
    for n in range(4, 7):
        for phase in ("solv", "gas"):
            elems, coords = read_xyz(BASE / f"mg_{lig}{n}_{phase}" / "final.xyz")
            mx, my, mz = coords[0]
            o_dists = [
                math.sqrt((cx - mx)**2 + (cy - my)**2 + (cz - mz)**2)
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
