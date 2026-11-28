import subprocess, re, json
from pathlib import Path

HARTREE_KCAL = 627.509474
BASE = Path("biohybrid_screen_v2")
OUT = Path("v3_thermo_hess")
OUT.mkdir(exist_ok=True)
Path("results").mkdir(exist_ok=True)

SPECIES = {
    "DME":             ("free_DME/xtbopt.xyz", 0),
    "Methoxyethanol":  ("free_Methoxyethanol/xtbopt.xyz", 0),
    "DMA":             ("free_Amide_DMA/xtbopt.xyz", 0),
    "Betaine":         ("free_Zwitterion_Bet/xtbopt.xyz", 0),
    "Mg_DME3":         ("complex_Mg_DME2_DME/xtbopt.xyz", 2),
    "Mg_DME2_MeOEt":   ("complex_Mg_DME2_Methoxyethanol/xtbopt.xyz", 2),
    "Mg_DME2_DMA":     ("complex_Mg_DME2_Amide_DMA/xtbopt.xyz", 2),
    "Mg_DME2_Betaine": ("complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2),
}

res_sp = {}
for name, (rel_p, chg) in SPECIES.items():
    w = OUT / name
    w.mkdir(exist_ok=True)
    (w / "geom.xyz").write_text((BASE / rel_p).read_text())
    p = subprocess.run(["xtb", "geom.xyz", "--ohess", "tight", "--alpb", "thf", "--chrg", str(chg)],
                       cwd=w, capture_output=True, text=True)
    txt = p.stdout + "\n" + p.stderr
    (w / "hess.out").write_text(txt)
    e = re.search(r"TOTAL ENERGY\s+([-\d.]+)\s+Eh", txt)
    g = re.search(r"TOTAL FREE ENERGY\s+([-\d.]+)\s+Eh", txt)
    z = re.search(r"zero point energy\s+([-\d.]+)\s+Eh", txt)
    res_sp[name] = {
        "E_Eh": float(e.group(1)) if e else None,
        "G_Eh": float(g.group(1)) if g else None,
        "ZPVE_Eh": float(z.group(1)) if z else None,
    }

rxns = {"Methoxyethanol": ("Mg_DME2_MeOEt", "Methoxyethanol"), "DMA": ("Mg_DME2_DMA", "DMA"), "Betaine": ("Mg_DME2_Betaine", "Betaine")}
exch = {}
for lig, (ck, lk) in rxns.items():
    dE = (res_sp[ck]["E_Eh"] + res_sp["DME"]["E_Eh"] - res_sp["Mg_DME3"]["E_Eh"] - res_sp[lk]["E_Eh"]) * HARTREE_KCAL
    dG = (res_sp[ck]["G_Eh"] + res_sp["DME"]["G_Eh"] - res_sp["Mg_DME3"]["G_Eh"] - res_sp[lk]["G_Eh"]) * HARTREE_KCAL
    exch[lig] = {"dE_kcal": round(dE, 2), "dG_kcal": round(dG, 2), "ddG_corr": round(dG - dE, 2)}

Path("results/v3_thermo.json").write_text(json.dumps({"species": res_sp, "exchange": exch}, indent=2))
print(json.dumps(exch, indent=2))
