import json
from pathlib import Path

# 1. Write clean v4_04_func_and_red.py
v4_04_code = '''import json
from pathlib import Path
from pyscf import dft, solvent
from pyscf.geomopt.geometric_solver import optimize
from v4_common import SPECIES, WORK, HARTREE_KCAL, HARTREE_EV, EPS_DME, load_mol, make_rks

opt = json.loads(Path("results/v4_dft_thermo.json").read_text())["species"]
funcs = [("b3lyp", "d3bj"), ("pbe0", "d3bj"), ("wb97x-d", None), ("m062x", None)]
enes = {f: {} for f, _ in funcs}

for name, item in SPECIES.items():
    _, chg = item
    mol = load_mol(opt[name]["xyz"], "def2-tzvp", chg, 0)
    for xc, disp in funcs:
        enes[xc][name] = float(make_rks(mol, xc, disp).kernel())

rxns = {"Methoxyethanol": ("Mg_DME2_MeOEt", "Methoxyethanol"), "DMA": ("Mg_DME2_DMA", "DMA"), "Betaine": ("Mg_DME2_Betaine", "Betaine")}
fsum = {lig: {xc: round((enes[xc][ck] + enes[xc]["DME"] - enes[xc]["Mg_DME3"] - enes[xc][lk]) * HARTREE_KCAL, 2)
              for xc, _ in funcs} for lig, (ck, lk) in rxns.items()}
Path("results/v4_functionals.json").write_text(json.dumps(fsum, indent=2))

def make_uks(mol):
    umf = dft.UKS(mol)
    umf.xc = "b3lyp"
    umf.disp = "d3bj"
    umf = solvent.PCM(umf)
    umf.with_solvent.method = "IEF-PCM"
    umf.with_solvent.eps = EPS_DME
    return umf

rdir = WORK / "04_red"
rdir.mkdir(exist_ok=True)
red_out = {}
for name in ["DME", "DMA", "Betaine", "Mg_DME3", "Mg_DME2_DMA", "Mg_DME2_Betaine"]:
    _, chg = SPECIES[name]
    e_ox = float(make_rks(load_mol(opt[name]["xyz"], "ma-def2-tzvp", chg, 0)).kernel())
    mol_ropt = optimize(make_uks(load_mol(opt[name]["xyz"], "def2-svp", chg - 1, 1)), maxsteps=60)
    rxyz = rdir / f"{name}_red.xyz"
    mol_ropt.tofile(str(rxyz), format="xyz")
    e_red = float(make_uks(load_mol(rxyz, "ma-def2-tzvp", chg - 1, 1)).kernel())
    red_out[name] = {"AEA_eV": round((e_ox - e_red) * HARTREE_EV, 3)}
    print(f"AEA {name}: {red_out[name]}")

Path("results/v4_adiabatic_red.json").write_text(json.dumps(red_out, indent=2))
'''
Path("v4_04_func_and_red.py").write_text(v4_04_code)

# 2. Write clean v4_05_neb.py using real methyl_C indices from Task 2
audit = json.loads(Path("results/v3_geom_audit.json").read_text())
pure_idx = audit.get("NEB_pure_start", {}).get("methyl_C_indices", [])
hyb_idx  = audit.get("NEB_hyb_start", {}).get("methyl_C_indices", [])

v4_05_code = f'''import json
from pathlib import Path
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.mep import NEB
from ase.optimize import FIRE
from xtb.ase.calculator import XTB

odir = Path("v4_work/05_neb")
odir.mkdir(exist_ok=True)

def run_neb(s_xyz, e_xyz, label, fixed_idx):
    ini, fin = read(s_xyz), read(e_xyz)
    imgs = [ini] + [ini.copy() for _ in range(7)] + [fin]
    for im in imgs:
        im.set_constraint(FixAtoms(indices=fixed_idx))
        im.calc = XTB(method="GFN2-xTB", solvent="THF", charge=2)
    neb = NEB(imgs, climb=True, k=0.15)
    neb.interpolate(method="idpp")
    FIRE(neb, trajectory=str(odir / f"{{label}}.traj")).run(fmax=0.05, steps=250)
    e0 = next(iter(imgs)).get_potential_energy()
    ec = [round((im.get_potential_energy() - e0) * 23.0605, 2) for im in imgs]
    write(str(odir / f"{{label}}_path.xyz"), imgs)
    return {{"label": label, "barrier_kcal": max(ec), "profile": ec}}

res = [
    run_neb("relaxed_hop_v2/pure_bet_fwd_0/xtbopt.xyz", "relaxed_hop_v2/pure_bet_fwd_14/xtbopt.xyz", "pure_bet", {pure_idx}),
    run_neb("relaxed_hop_v2/hybrid_dma_fwd_0/xtbopt.xyz", "relaxed_hop_v2/hybrid_dma_fwd_14/xtbopt.xyz", "hybrid_dma", {hyb_idx}),
]
Path("results/v4_neb.json").write_text(json.dumps(res, indent=2))
'''
Path("v4_05_neb.py").write_text(v4_05_code)

# 3. Write clean run_v4_all.sh
run_sh = """#!/usr/bin/env bash
set -e
ulimit -s unlimited
bash v4_01_setup.sh
conda run -n mg_nanobio python3 v4_02_crest.py
conda run -n mg_nanobio python3 v4_03_dft_opt.py
conda run -n mg_nanobio python3 v4_04_func_and_red.py
conda run -n mg_nanobio python3 v4_05_neb.py
echo "ALL_V4_COMPLETE"
"""
Path("run_v4_all.sh").write_text(run_sh)
print("PATCHED_OK")
