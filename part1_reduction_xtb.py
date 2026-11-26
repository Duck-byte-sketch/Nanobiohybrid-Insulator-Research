import json, shutil, subprocess
from pathlib import Path

HARTREE_TO_KCAL = 627.509474
HARTREE_TO_EV = 27.211386245988
MG_REF_SHIFT_V = 2.07

W = Path("biohybrid_screen_v2")
W_RED = Path("reduction_alpb_thf")
W_RED.mkdir(exist_ok=True)
RES_DIR = Path("results")
RES_DIR.mkdir(exist_ok=True)
XTB_BIN = shutil.which("xtb")

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

def run_xtb_sp(xyz_path, charge, uhf, tag):
    d = W_RED / tag
    d.mkdir(exist_ok=True)
    shutil.copy(xyz_path, d / "geom.xyz")
    cmd = [XTB_BIN, "geom.xyz", "--sp", "--chrg", str(charge), "--uhf", str(uhf), "--alpb", "thf"]
    res = subprocess.run(cmd, cwd=d, capture_output=True, text=True)
    out_all = res.stdout + "\n" + res.stderr
    (d / "xtb.out").write_text(out_all)
    e = get_energy(out_all)
    spin_info = {}
    if uhf == 1:
        atoms = []
        for idx, l in enumerate(Path(xyz_path).read_text().splitlines()):
            if idx >= 2 and l.strip():
                sym, *_ = l.split()
                atoms.append(sym)
        sp_file = d / "spinpop"
        if sp_file.exists():
            sp_vals = [float(x.split().pop()) for x in sp_file.read_text().splitlines() if x.strip()]
            spins = list(zip(atoms, sp_vals))
            mg_sp = sum(v for s, v in spins if s == "Mg")
            lig_sp = sum(v for s, v in spins if s != "Mg")
            top = sorted(spins, key=lambda item: abs(item.pop()) if False else abs(list(item).pop()), reverse=True)[:4]
            spin_info = {
                "mg_spin": round(mg_sp, 3),
                "ligand_spin": round(lig_sp, 3),
                "top_atoms": [(s, round(v, 3)) for s, v in top],
            }
    return e, spin_info

print("=== PART 1: VERTICAL ELECTRON REDUCTION & SPIN LOCALIZATION (GFN2-xTB / ALPB-THF) ===", flush=True)
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
    safe = name.replace("[", "").replace("]", "").replace("+", "p")
    e_ox, _ = run_xtb_sp(W / rel_xyz, q0, 0, f"{safe}_ox")
    e_red, sp_inf = run_xtb_sp(W / rel_xyz, q0 - 1, 1, f"{safe}_red")
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
for name, d in reduction_results.items():
    shift = d["vea_ev"] - ctrl_vea if d["charge_ox"] == 2 else None
    d["shift_vs_MgDME3_eV"] = round(shift, 3) if shift is not None else None
    print(f"{name:27s} | VEA={d['vea_ev']:+6.2f} eV | E_red(vs Mg/Mg2+)~{d['approx_E_red_vs_Mg_V']:+6.2f} V | shift={str(d['shift_vs_MgDME3_eV']):>6s} | spin={d['spin_localization']}", flush=True)

payload = {
    "swap_energies": SWAP_SUMMARY,
    "tzvp_pcm_raw_eh": TZVP_PCM_ENERGIES,
    "vertical_reduction_alpb_thf": reduction_results,
}
(RES_DIR / "verified_v2_data.json").write_text(json.dumps(payload, indent=2) + "\n")
print("Saved Part 1 cleanly to results/verified_v2_data.json!", flush=True)
