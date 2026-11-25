from pathlib import Path
from pyscf import gto, dft

HARTREE_TO_KCAL = 627.509474
W = Path("biohybrid_screen")

def xyz_to_atoms(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines.pop(0).strip())
    lines.pop(0)
    return "\n".join(lines[:n])

def run_pyscf_direct(xyz_path, charge, label):
    mol = gto.M(atom=xyz_to_atoms(xyz_path), basis="def2-svp", charge=charge, spin=0, verbose=0)
    mf = dft.RKS(mol)
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf.max_memory = 700
    mf.direct_scf = True
    mf.direct_scf_tol = 1e-10
    mf.grids.level = 2
    mf.conv_tol = 1e-6
    e = mf.kernel()
    print(f"[DFT B3LYP-D3BJ/def2-SVP] {label:22s} : {e:16.8f} Eh | converged={mf.converged}", flush=True)
    return e

print("=== PHASE 3: DFT VALIDATION OF ZWITTERIONIC BETAINE DME-DISPLACEMENT ===", flush=True)
e_dme = run_pyscf_direct(W / "free_DME/xtbopt.xyz", 0, "DME")
e_bet = run_pyscf_direct(W / "free_Zwitterion_Bet/xtbopt.xyz", 0, "Zwitterion_Bet")
e_mg_dme3 = run_pyscf_direct(W / "complex_Mg_DME2_DME/xtbopt.xyz", 2, "[Mg(DME)3]2+")
e_mg_dme2_bet = run_pyscf_direct(W / "complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2, "[Mg(DME)2(Bet)]2+")

de_swap_dft = ((e_mg_dme2_bet + e_dme) - (e_mg_dme3 + e_bet)) * HARTREE_TO_KCAL
print(f"\n[RESULT] Gas-Phase DFT dE_swap ([Mg(DME)3]2+ + Betaine -> [Mg(DME)2(Bet)]2+ + DME): {de_swap_dft:8.2f} kcal/mol", flush=True)
