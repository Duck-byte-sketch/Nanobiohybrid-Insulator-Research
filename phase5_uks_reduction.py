import json
from pathlib import Path
import numpy as np
from pyscf import dft, gto, solvent

HARTREE_TO_EV = 27.211386245988
W_V2 = Path("biohybrid_screen_v2")
RES_DIR = Path("results")
RES_DIR.mkdir(exist_ok=True)


def load_xyz_str(path):
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    it = iter(lines)
    n = int(next(it).strip())
    next(it)
    body = []
    for _ in range(n):
        body.append(next(it))
    return "\n".join(body)


def run_single_point(name, xyz_str, charge, spin, basis="ma-def2-tzvp", eps=7.58):
    mol = gto.M(
        atom=xyz_str,
        basis=basis,
        charge=charge,
        spin=spin,
        verbose=0,
    )
    mf = dft.RKS(mol) if spin == 0 else dft.UKS(mol)
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf = mf.density_fit(auxbasis="def2-universal-jkfit")
    mf = solvent.PCM(mf)
    mf.with_solvent.method = "IEF-PCM"
    mf.with_solvent.eps = eps
    mf.conv_tol = 1e-7
    mf.max_cycle = 150
    e_tot = float(mf.kernel())

    spin_info = {}
    if spin > 0:
        s2, mult = mf.spin_square()
        dm_alpha, dm_beta = mf.make_rdm1()
        p_spin = dm_alpha - dm_beta
        s_ovlp = mol.intor("int1e_ovlp")
        atom_spins = []
        aoslices = mol.aoslice_by_atom()
        for ia in range(mol.natm):
            _, _, p0, p1 = aoslices[ia]
            sp_val = float(np.einsum("ij,ji->", p_spin[p0:p1, :], s_ovlp[:, p0:p1]))
            atom_spins.append((mol.atom_symbol(ia), round(sp_val, 4)))
        mg_spin = round(sum(v for s, v in atom_spins if s == "Mg"), 4)
        lig_spin = round(sum(v for s, v in atom_spins if s != "Mg"), 4)
        top_atoms = sorted(atom_spins, key=lambda x: abs(x[1]), reverse=True)[:5]
        spin_info = {
            "S2": round(float(s2), 4),
            "multiplicity": round(float(mult), 4),
            "mg_spin_mulliken": mg_spin,
            "ligand_spin_mulliken": lig_spin,
            "total_spin_check": round(mg_spin + lig_spin, 4),
            "top_spin_atoms": top_atoms,
        }
    return {
        "converged": bool(mf.converged),
        "energy_eh": round(e_tot, 8),
        "spin_analysis": spin_info,
    }


def main():
    targets = (
        ("Free_DME_Calib", W_V2 / "free_DME/xtbopt.xyz", 0),
        ("Free_DMA_Calib", W_V2 / "free_Amide_DMA/xtbopt.xyz", 0),
        ("Free_Betaine_Calib", W_V2 / "free_Zwitterion_Bet/xtbopt.xyz", 0),
        ("Mg(DME)3_2+", W_V2 / "Mg_DME3/xtbopt.xyz", 2),
        ("Mg(DME)2(DMA)_2+", W_V2 / "Mg_DME2_Amide_DMA/xtbopt.xyz", 2),
        ("Mg(DME)2(Betaine)_2+", W_V2 / "Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2),
    )

    ec_xyz = """C   0.000000   0.000000   0.768500
O   0.000000   0.000000   1.976200
O   0.000000   1.121500  -0.011200
O   0.000000  -1.121500  -0.011200
C   0.000000   0.765400  -1.391000
C   0.000000  -0.765400  -1.391000
H   0.892100   1.201400  -1.852000
H  -0.892100   1.201400  -1.852000
H   0.892100  -1.201400  -1.852000
H  -0.892100  -1.201400  -1.852000"""

    results = {}
    print("=== UKS B3LYP-D3BJ/ma-def2-TZVP/IEF-PCM(THF) VERTICAL REDUCTION ===", flush=True)

    all_jobs = [("EC_Ref_Calib", ec_xyz, 0)]
    for name, path, q in targets:
        if path.exists():
            all_jobs.append((name, load_xyz_str(path), q))

    for name, xyz_str, q_ox in all_jobs:
        ox = run_single_point(name + "_ox", xyz_str, charge=q_ox, spin=0)
        red = run_single_point(name + "_red", xyz_str, charge=q_ox - 1, spin=1)
        vea_eh = ox["energy_eh"] - red["energy_eh"]
        vea_ev = round(vea_eh * HARTREE_TO_EV, 3)
        e_vs_mg = round(vea_ev - 2.07, 3)

        sp = red["spin_analysis"]
        results[name] = {
            "charge_ox": q_ox,
            "charge_red": q_ox - 1,
            "E_ox_eh": ox["energy_eh"],
            "E_red_eh": red["energy_eh"],
            "VEA_eV": vea_ev,
            "E_red_est_V_vs_Mg": e_vs_mg,
            "spin_localization": sp,
        }
        print(
            f"{name:22s} | VEA={vea_ev:+6.3f} eV | E_vs_Mg={e_vs_mg:+6.3f} V | "
            f"Mg_spin={sp.get('mg_spin_mulliken', 0.0):+6.3f} | Lig_spin={sp.get('ligand_spin_mulliken', 0.0):+6.3f}",
            flush=True,
        )

    out_path = RES_DIR / "uks_reduction_benchmark.json"
    out_path.write_text(json.dumps({
        "method": "RKS/UKS B3LYP-D3BJ / ma-def2-TZVP / IEF-PCM(THF, eps=7.58)",
        "reference_conversion": "E_red(V vs Mg2+/Mg) ~ VEA(eV) - 2.07 V",
        "species": results,
    }, indent=2), encoding="utf-8")

    v2_path = RES_DIR / "verified_v2_data.json"
    if v2_path.exists():
        v2 = json.loads(v2_path.read_text(encoding="utf-8"))
        v2["uks_b3lyp_ma_tzvp_reduction"] = results
        v2_path.write_text(json.dumps(v2, indent=2), encoding="utf-8")
    print(f"\nSaved UKS reduction benchmark to {out_path}")


if __name__ == "__main__":
    main()
