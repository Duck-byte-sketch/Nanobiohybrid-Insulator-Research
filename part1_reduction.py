import gc, json, os
from pathlib import Path
import numpy as np
from pyscf import gto, dft

os.environ["PYSCF_TMPDIR"] = os.path.expanduser("~/pyscf_scratch")
Path(os.environ["PYSCF_TMPDIR"]).mkdir(parents=True, exist_ok=True)

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_EV = 27.211386245988
MG_REF_SHIFT_V = 2.07
W = Path("biohybrid_screen_v2")
RES_DIR = Path("results")
RES_DIR.mkdir(exist_ok=True)

TZVP_PCM_ENERGIES = {
    "Free_DME": -309.00994568,
    "Free_Methoxyethanol": -269.69373563,
    "Free_Amide_DMA": -287.97617863,
    "Free_Zwitterion_Bet": -402.54246814,
    "Complex_Mg_DME3": -1127.00960969,
    "Complex_Mg_DME2_MeOEtOH": -1087.66309524,
    "Complex_Mg_DME2_DMA": -1105.96045504,
    "Complex_Mg_DME2_Bet": -1220.55739072,
}

SWAP_SUMMARY = {
    "Methoxyethanol": {
        "xtb_alpb_kcal": 23.41,
        "dft_tzvp_pcm_kcal": round(((TZVP_PCM_ENERGIES["Complex_Mg_DME2_MeOEtOH"] + TZVP_PCM_ENERGIES["Free_DME"]) -
                                    (TZVP_PCM_ENERGIES["Complex_Mg_DME3"] + TZVP_PCM_ENERGIES["Free_Methoxyethanol"])) * HARTREE_TO_KCAL, 2),
        "cn": 5,
    },
    "DMA_amide": {
        "xtb_alpb_kcal": 10.77,
        "dft_tzvp_pcm_kcal": round(((TZVP_PCM_ENERGIES["Complex_Mg_DME2_DMA"] + TZVP_PCM_ENERGIES["Free_DME"]) -
                                    (TZVP_PCM_ENERGIES["Complex_Mg_DME3"] + TZVP_PCM_ENERGIES["Free_Amide_DMA"])) * HARTREE_TO_KCAL, 2),
        "cn": 5,
    },
    "Betaine": {
        "xtb_alpb_kcal": -12.35,
        "dft_tzvp_pcm_kcal": round(((TZVP_PCM_ENERGIES["Complex_Mg_DME2_Bet"] + TZVP_PCM_ENERGIES["Free_DME"]) -
                                    (TZVP_PCM_ENERGIES["Complex_Mg_DME3"] + TZVP_PCM_ENERGIES["Free_Zwitterion_Bet"])) * HARTREE_TO_KCAL, 2),
        "cn": 6,
    },
}

def xyz_to_atoms_str(xyz_path):
    lines = Path(xyz_path).read_text().strip().splitlines()
    n = int(lines[0].strip())
    return "\n".join(lines[2:2+n])

def run_svp_pcm_sp(xyz_path, charge, spin, label):
    mol = gto.M(atom=xyz_to_atoms_str(xyz_path), basis="def2-svp", charge=charge, spin=spin, verbose=0)
    mf = (dft.RKS(mol) if spin == 0 else dft.UKS(mol)).density_fit(auxbasis="def2-universal-jkfit")
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf.max_memory = 3500
    mf.conv_tol = 1e-6
    mf = mf.PCM()
    mf.with_solvent.method = "IEF-PCM"
    mf.with_solvent.eps = 7.58
    e = float(mf.kernel())
    conv = bool(mf.converged)

    spin_info = {}
    if spin == 1:
        dm_a, dm_b = mf.make_rdm1()
        s_ovlp = mol.intor("int1e_ovlp")
        pop_a = np.einsum("ij,ji->i", dm_a, s_ovlp)
        pop_b = np.einsum("ij,ji->i", dm_b, s_ovlp)
        ao_slices = mol.aoslice_by_atom()
        atom_spins = []
        for ia in range(mol.natm):
            p0, p1 = ao_slices[ia, 2], ao_slices[ia, 3]
            sp_val = float(np.sum(pop_a[p0:p1] - pop_b[p0:p1]))
            atom_spins.append((mol.atom_symbol(ia), sp_val))
        mg_spin = sum(v for sym, v in atom_spins if sym == "Mg")
        lig_spin = sum(v for sym, v in atom_spins if sym != "Mg")
        top_atoms = sorted(atom_spins, key=lambda x: abs(x), reverse=True)[:4]
        spin_info = {
            "mg_spin": round(mg_spin, 3),
            "ligand_spin": round(lig_spin, 3),
            "top_atoms": [(sym, round(val, 3)) for sym, val in top_atoms],
        }

    del mf, mol
    gc.collect()
    print(f"[SVP-PCM(THF)] {label:26s} (q={charge:+d}, S={spin}) : {e:15.8f} Eh | conv={conv} | spin={spin_info}", flush=True)
    return e, conv, spin_info

print("=== PART 1: VERTICAL ELECTRON REDUCTION & SPIN LOCALIZATION (IEF-PCM THF) ===", flush=True)
TARGETS = [
    ("Control_[Mg(DME)3]2+",      "complex_Mg_DME2_DME/xtbopt.xyz",            2),
    ("Complex_[Mg(DME)2(DMA)]2+", "complex_Mg_DME2_Amide_DMA/xtbopt.xyz",      2),
    ("Complex_[Mg(DME)2(Bet)]2+", "complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2),
    ("Neutral_Free_DME",          "free_DME/xtbopt.xyz",                       0),
    ("Neutral_Free_DMA",          "free_Amide_DMA/xtbopt.xyz",                 0),
    ("Neutral_Free_Betaine",      "free_Zwitterion_Bet/xtbopt.xyz",            0),
]

reduction_results = {}
for name, rel_xyz, q0 in TARGETS:
    e_ox, _, _ = run_svp_pcm_sp(W / rel_xyz, q0, 0, f"{name} (ox)")
    e_red, _, sp_inf = run_svp_pcm_sp(W / rel_xyz, q0 - 1, 1, f"{name} (red)")
    vea_ev = (e_ox - e_red) * HARTREE_TO_EV
    ered_mg = vea_ev - MG_REF_SHIFT_V
    reduction_results[name] = {
        "charge_ox": q0,
        "e_ox_eh": round(e_ox, 8),
        "e_red_eh": round(e_red, 8),
        "vea_ev": round(vea_ev, 3),
        "approx_E_red_vs_Mg_V": round(ered_mg, 3),
        "spin_localization": sp_inf,
    }

ctrl_vea = reduction_results["Control_[Mg(DME)3]2+"]["vea_ev"]
print("\n--- Summary vs. [Mg(DME)3]2+ Baseline Control ---", flush=True)
for name, d in reduction_results.items():
    shift = d["vea_ev"] - ctrl_vea if d["charge_ox"] == 2 else None
    d["shift_vs_MgDME3_eV"] = round(shift, 3) if shift is not None else None
    print(f"{name:27s} | VEA={d['vea_ev']:+6.2f} eV | E_red(vs Mg/Mg2+)~{d['approx_E_red_vs_Mg_V']:+6.2f} V | shift={str(d['shift_vs_MgDME3_eV']):>6s} | spin={d['spin_localization']}", flush=True)

payload = {
    "swap_energies": SWAP_SUMMARY,
    "tzvp_pcm_raw_eh": TZVP_PCM_ENERGIES,
    "vertical_reduction_svp_pcm": reduction_results,
}
(RES_DIR / "verified_v2_data.json").write_text(json.dumps(payload, indent=2) + "\n")
print("Saved Part 1 to results/verified_v2_data.json!", flush=True)
