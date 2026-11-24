import math
from pathlib import Path

HARTREE_TO_KCAL = 627.509474

def read_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines[0].strip())
    elems, coords = [], []
    for line in lines[2:2+n]:
        parts = line.split()
        elems.append(parts[0])
        coords.append((float(parts[1]), float(parts[2]), float(parts[3])))
    return elems, coords

def get_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            return float(line.split()[3])
    raise ValueError(f"No TOTAL ENERGY in {out_path}")

base = Path("mg_benchmark")
print("=== 1. Coordination Audit (Mg-O < 2.4 A) ===")
for lig in ("dme", "urea"):
    for n in range(4, 7):
        for phase in ("solv", "gas"):
            xyz = base / f"mg_{lig}{n}_{phase}" / "final.xyz"
            elems, coords = read_xyz(xyz)
            mx, my, mz = coords[0]
            o_dists = [
                math.sqrt((cx - mx)**2 + (cy - my)**2 + (cz - mz)**2)
                for e, (cx, cy, cz) in zip(elems, coords) if e == "O"
            ]
            cn = sum(d < 2.40 for d in o_dists)
            near = sorted(round(d, 2) for d in o_dists)[:8]
            print(f"mg_{lig}{n}_{phase:4s} | CN(Mg-O)={cn} | nearest O (A): {near}")

print("\n=== 2. Stepwise Desolvation Energies (kcal/mol) ===")
for phase in ("solv", "gas"):
    print(f"\n--- Phase: {phase.upper()} ---")
    for lig in ("dme", "urea"):
        e_lig = get_energy(base / f"{lig}_{phase}" / "final.out")
        e4 = get_energy(base / f"mg_{lig}4_{phase}" / "final.out")
        e5 = get_energy(base / f"mg_{lig}5_{phase}" / "final.out")
        e6 = get_energy(base / f"mg_{lig}6_{phase}" / "final.out")
        de_6_5 = ((e5 + e_lig) - e6) * HARTREE_TO_KCAL
        de_5_4 = ((e4 + e_lig) - e5) * HARTREE_TO_KCAL
        print(f"{lig.upper():4s} | n=6->5: {de_6_5:8.2f} kcal/mol | n=5->4: {de_5_4:8.2f} kcal/mol")
