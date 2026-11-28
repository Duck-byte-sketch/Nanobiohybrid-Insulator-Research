import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"

import sys
import json
import time
from pathlib import Path
import numpy as np
from pyscf import dft, gto

HARTREE_TO_EV = 27.211386245988
RES_DIR = Path("results")


def load_xyz_str(path):
    lines = Path(path).read_text(encoding="utf-8").strip().splitlines()
    it = iter(lines)
    n = int(next(it).strip())
    next(it)
    return "\n".join(next(it) for _ in range(n))


def run_lowmem_sp(xyz_str, charge, spin, basis="def2-svp", eps=7.43):
    mol = gto.M(atom=xyz_str, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = dft.RKS(mol) if spin == 0 else dft.UKS(mol)
    mf.max_memory = 1500
    mf.direct_scf = True
    mf.xc = "b3lyp"
    mf.disp = "d3bj"
    mf = mf.ddCOSMO()
    mf.with_solvent.eps = eps
    mf.conv_tol = 1e-6
    mf.max_cycle = 120
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
        top_atoms = sorted(atom_spins, key=lambda x: abs(list(x).pop()), reverse=True)[:5]
        spin_info = {
            "S2": round(float(s2), 4),
            "multiplicity": round(float(mult), 4),
            "mg_spin_mulliken": mg_spin,
            "ligand_spin_mulliken": lig_spin,
            "top_spin_atoms": top_atoms,
        }
    return {"converged": bool(mf.converged), "energy_eh": round(e_tot, 8), "spin_analysis": spin_info}


def main():
    it_arg = iter(sys.argv)
    next(it_arg)
    target_filter = next(it_arg, "")

    jobs = (
        ("Mg_DME3", "Mg(DME)3_2+", "biohybrid_screen_v2/complex_Mg_DME2_DME/xtbopt.xyz", 2),
        ("Mg_DMA", "Mg(DME)2(DMA)_2+", "biohybrid_screen_v2/complex_Mg_DME2_Amide_DMA/xtbopt.xyz", 2),
        ("Mg_Bet", "Mg(DME)2(Betaine)_2+", "biohybrid_screen_v2/complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz", 2),
    )

    out_path = RES_DIR / "uks_reduction_benchmark.json"
    bench = json.loads(out_path.read_text(encoding="utf-8"))
    results = bench.get("species", {})

    for key, label, xyz_rel, q_ox in jobs:
        if target_filter and target_filter != key:
            continue
        t0 = time.time()
        xyz_str = load_xyz_str(xyz_rel)
        ox = run_lowmem_sp(xyz_str, charge=q_ox, spin=0)
        red = run_lowmem_sp(xyz_str, charge=q_ox - 1, spin=1)
        vea_ev = round((ox["energy_eh"] - red["energy_eh"]) * HARTREE_TO_EV, 3)
        e_vs_mg = round(vea_ev - 2.07, 3)
        sp = red["spin_analysis"]
        dt = round(time.time() - t0, 1)

        results[label] = {
            "basis": "def2-SVP",
            "solvation": "COSMO(eps=7.43, lebedev=11)",
            "charge_ox": q_ox,
            "charge_red": q_ox - 1,
            "E_ox_eh": ox["energy_eh"],
            "E_red_eh": red["energy_eh"],
            "VEA_eV": vea_ev,
            "E_red_est_V_vs_Mg": e_vs_mg,
            "spin_localization": sp,
        }
        bench["species"] = results
        out_path.write_text(json.dumps(bench, indent=2), encoding="utf-8")

        v2_path = RES_DIR / "verified_v2_data.json"
        if v2_path.exists():
            v2 = json.loads(v2_path.read_text(encoding="utf-8"))
            v2["uks_b3lyp_reduction"] = results
            v2_path.write_text(json.dumps(v2, indent=2), encoding="utf-8")

        print(
            f"{label:22s} ({dt:5.1f}s) | VEA={vea_ev:+6.3f} eV | E_vs_Mg={e_vs_mg:+6.3f} V | "
            f"Mg_spin={sp.get('mg_spin_mulliken', 0.0):+6.3f} | Lig_spin={sp.get('ligand_spin_mulliken', 0.0):+6.3f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
