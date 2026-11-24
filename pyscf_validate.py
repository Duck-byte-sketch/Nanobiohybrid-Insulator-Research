from pathlib import Path
from pyscf import gto, dft

HARTREE_TO_KCAL = 627.509474
BASE = Path("mg_benchmark")

def get_xtb_energy(out_path):
    for line in reversed(Path(out_path).read_text().splitlines()):
        if "TOTAL ENERGY" in line:
            for tok in line.split():
                try:
                    val = float(tok)
                    if val < -1.0:
                        return val
                except ValueError:
                    continue
    raise ValueError(f"Could not parse TOTAL ENERGY in {out_path}")

def run_pyscf_direct(atom_spec, charge, label):
    mol = gto.M(atom=atom_spec, basis="def2-svp", charge=charge, spin=0, verbose=0)
    # Pure Direct SCF (no density_fit -> zero 3 GiB array, <150 MB RAM)
    mf = dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf.max_memory = 700
    mf.direct_scf = True
    mf.direct_scf_tol = 1e-10
    mf.grids.level = 2
    mf.conv_tol = 1e-6
    e_dft = mf.kernel()
    print(f"[DFT B3LYP-D3BJ/def2-SVP] {label:16s} : {e_dft:16.8f} Eh | converged={mf.converged}", flush=True)
    return e_dft

def xyz_to_atoms(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    return "\n".join(lines[:n])

# Already converged in previous step:
e_mg_dft   = -199.17113857
e_dme_dft  = -308.64183326
e_urea_dft = -225.11031106

print("=== PYSCF DFT SINGLE-POINT VALIDATION (DIRECT SCF, GAS PHASE, n=4) ===", flush=True)
e_mg_urea4_dft = run_pyscf_direct(xyz_to_atoms(BASE / "mg_urea4_gas/final.xyz"), 2, "[Mg(Urea)4]2+")
e_mg_dme4_dft  = run_pyscf_direct(xyz_to_atoms(BASE / "mg_dme4_gas/final.xyz"), 2, "[Mg(DME)4]2+")

bind_dme4_dft  = (e_mg_dme4_dft  - (e_mg_dft + 4 * e_dme_dft))  * HARTREE_TO_KCAL
bind_urea4_dft = (e_mg_urea4_dft - (e_mg_dft + 4 * e_urea_dft)) * HARTREE_TO_KCAL

e_dme_xtb      = get_xtb_energy(BASE / "dme_gas/final.out")
e_urea_xtb     = get_xtb_energy(BASE / "urea_gas/final.out")
e_mg_dme4_xtb  = get_xtb_energy(BASE / "mg_dme4_gas/final.out")
e_mg_urea4_xtb = get_xtb_energy(BASE / "mg_urea4_gas/final.out")

exchange_dft = ((e_mg_urea4_dft + 4 * e_dme_dft)  - (e_mg_dme4_dft + 4 * e_urea_dft))  * HARTREE_TO_KCAL
exchange_xtb = ((e_mg_urea4_xtb + 4 * e_dme_xtb) - (e_mg_dme4_xtb + 4 * e_urea_xtb)) * HARTREE_TO_KCAL

print("\n=== BENCHMARK SUMMARY: GFN2-xTB vs. B3LYP-D3BJ/def2-SVP ===", flush=True)
print(f"DFT Total Formation Energy [Mg(DME)4]2+  : {bind_dme4_dft:8.2f} kcal/mol")
print(f"DFT Total Formation Energy [Mg(Urea)4]2+ : {bind_urea4_dft:8.2f} kcal/mol")
print(f"Ligand Exchange ([Mg(DME)4]2+ + 4 Urea -> [Mg(Urea)4]2+ + 4 DME):")
print(f"  GFN2-xTB : {exchange_xtb:8.2f} kcal/mol")
print(f"  DFT      : {exchange_dft:8.2f} kcal/mol")
