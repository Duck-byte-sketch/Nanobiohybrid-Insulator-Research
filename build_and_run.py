import math, subprocess, shutil, sys
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem

HARTREE_TO_KCAL = 627.509474
BASE = Path("mg_benchmark")

xtb_bin = shutil.which("xtb")
if not xtb_bin:
    sys.exit("ERROR: xtb not found on PATH.")
print(f"Using xtb binary: {xtb_bin}")

def write_xyz(path, elems, coords, comment=""):
    with open(path, "w") as f:
        f.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

def build_chelated_mg_dme3():
    smi = "CO1->[Mg+2]23(<-O(C)CC1)(<-O(C)CCO2C)<-O(C)CCO3C"
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        sys.exit("ERROR: RDKit failed to parse dative chelate SMILES.")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    AllChem.EmbedMolecule(mol, params)
    AllChem.UFFOptimizeMolecule(mol, maxIters=1000)

    conf = mol.GetConformer()
    mg_idx = next(a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "Mg")
    mp = conf.GetAtomPosition(mg_idx)
    mx, my, mz = mp.x, mp.y, mp.z

    rw = Chem.RWMol(mol)
    rw.RemoveAtom(mg_idx)
    frags = Chem.GetMolFrags(rw, asMols=True, sanitizeFrags=False)

    elems = ["Mg"]
    coords = [(0.0, 0.0, 0.0)]
    for frag in frags:
        fconf = frag.GetConformer()
        for atom in frag.GetAtoms():
            p = fconf.GetAtomPosition(atom.GetIdx())
            elems.append(atom.GetSymbol())
            coords.append((p.x - mx, p.y - my, p.z - mz))
    return elems, coords

def build_free_dme():
    mol = Chem.AddHs(Chem.MolFromSmiles("COCCOC"))
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    cx = sum(p.x for p in pts) / len(pts)
    cy = sum(p.y for p in pts) / len(pts)
    cz = sum(p.z for p in pts) / len(pts)
    elems = [a.GetSymbol() for a in mol.GetAtoms()]
    coords = [(p.x - cx, p.y - cy, p.z - cz) for p in pts]
    return elems, coords

core_e, core_c = build_chelated_mg_dme3()
dme1_e, dme1_c = build_free_dme()

outer_shifts = [
    ( 5.6,  0.0,  0.0),
    (-5.6,  0.0,  0.0),
    ( 0.0,  5.6,  0.0),
]

for n in range(4, 7):
    elems = list(core_e)
    coords = list(core_c)
    for k in range(n - 3):
        sx, sy, sz = outer_shifts[k]
        elems.extend(dme1_e)
        coords.extend([(x + sx, y + sy, z + sz) for x, y, z in dme1_c])
    write_xyz(BASE / f"mg_dme{n}_chelate.xyz", elems, coords, f"Tris-chelated [Mg(DME)_{n}]2+")

def run_xtb(work_dir, in_file, chrg, solv_flag):
    shutil.copy(in_file, work_dir / "input.xyz")
    cmd1 = [xtb_bin, "input.xyz", "--opt", "normal", "--gfn", "2", "--chrg", str(chrg), "--etemp", "500"] + solv_flag
    res1 = subprocess.run(cmd1, cwd=work_dir, capture_output=True, text=True)
    if res1.returncode != 0 or not (work_dir / "xtbopt.xyz").exists():
        log = res1.stdout + "\n" + res1.stderr
        print(f"[ERROR] Stage 1 xtb failed in {work_dir}:\n" + "\n".join(log.splitlines()[-20:]))
        sys.exit(1)

    cmd2 = [xtb_bin, "xtbopt.xyz", "--opt", "tight", "--gfn", "2", "--chrg", str(chrg)] + solv_flag
    res2 = subprocess.run(cmd2, cwd=work_dir, capture_output=True, text=True)
    if res2.returncode != 0 or not (work_dir / "xtbopt.xyz").exists():
        (work_dir / "final.out").write_text(res1.stdout + "\n" + res1.stderr)
    else:
        (work_dir / "final.out").write_text(res2.stdout + "\n" + res2.stderr)
    shutil.copy(work_dir / "xtbopt.xyz", work_dir / "final.xyz")

print("Running GFN2-xTB optimizations on DME clusters (n=4, 5, 6)...")
for n in range(4, 7):
    run_xtb(BASE / f"mg_dme{n}_solv", BASE / f"mg_dme{n}_chelate.xyz", 2, ["--alpb", "thf"])
    run_xtb(BASE / f"mg_dme{n}_gas",  BASE / f"mg_dme{n}_solv/final.xyz", 2, [])

def read_xyz(path):
    lines = Path(path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    elems, coords = [], []
    for line in lines[:n]:
        el, xs, ys, zs = line.split()
        elems.append(el)
        coords.append((float(xs), float(ys), float(zs)))
    return elems, coords

def get_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            _, _, _, val, _ = line.split()
            return float(val)
    raise ValueError(f"No TOTAL ENERGY in {out_path}")

print("\n=== 1. Verified Coordination Numbers (Mg-O < 2.4 A) ===")
for lig in ("dme", "urea"):
    for n in range(4, 7):
        for phase in ("solv", "gas"):
            elems, coords = read_xyz(BASE / f"mg_{lig}{n}_{phase}" / "final.xyz")
            mx, my, mz = coords.pop(0)
            elems.pop(0)
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
