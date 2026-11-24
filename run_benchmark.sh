#!/usr/bin/env bash
set -eo pipefail

if command -v conda >/dev/null 2>&1; then
    CONDA_BASE=$(conda info --base)
    source "$CONDA_BASE/etc/profile.d/conda.sh"
    conda activate mg_bench
else
    for p in "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "/opt/conda"; do
        if [ -f "$p/etc/profile.d/conda.sh" ]; then
            source "$p/etc/profile.d/conda.sh"
            conda activate mg_bench
            break
        fi
    done
fi

ulimit -s unlimited 2>/dev/null || true
export OMP_STACKSIZE=2G
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

WORKDIR="$(pwd)/mg_benchmark"
mkdir -p "$WORKDIR"
cd "$WORKDIR"

python3 - << 'PYEOF'
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

def make_mol_xyz(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    coords -= coords.mean(axis=0)
    elems = [a.GetSymbol() for a in mol.GetAtoms()]
    return elems, coords

def write_xyz(filename, elems, coords, comment=""):
    with open(filename, "w") as f:
        f.write(f"{len(elems)}\n{comment}\n")
        for e, (x, y, z) in zip(elems, coords):
            f.write(f"{e:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")

dme_e, dme_c = make_mol_xyz("COCCOC")
urea_e, urea_c = make_mol_xyz("NC(=O)N")
write_xyz("dme_init.xyz", dme_e, dme_c, "DME")
write_xyz("urea_init.xyz", urea_e, urea_c, "Urea")

oct_dirs = np.array([
    ( 1.0,  0.0,  0.0),
    (-1.0,  0.0,  0.0),
    ( 0.0,  1.0,  0.0),
    ( 0.0, -1.0,  0.0),
    ( 0.0,  0.0,  1.0),
    ( 0.0,  0.0, -1.0)
], dtype=float)

for name, (elems, coords), max_first in (("dme", (dme_e, dme_c), 3), ("urea", (urea_e, urea_c), 6)):
    for n in range(4, 7):
        all_e = ["Mg"]
        all_c = [np.zeros(3)]
        for i in range(n):
            shell_r = 3.2 if i < max_first else 6.2
            offset = np.array((0.5, 0.5, 0.5)) if i >= 6 else 0.0
            vec = oct_dirs[i % 6] * shell_r + offset
            all_e.extend(elems)
            all_c.extend(coords + vec)
        write_xyz(f"mg_{name}{n}_init.xyz", all_e, np.array(all_c), f"[Mg({name})_{n}]2+")

print("Initial structures written successfully.")
PYEOF

run_xtb_opt() {
    local name="$1"
    local xyz_in="$2"
    local chrg="$3"
    local solv_flag="$4"

    mkdir -p "$WORKDIR/$name"
    cp "$xyz_in" "$WORKDIR/$name/input.xyz"
    pushd "$WORKDIR/$name" > /dev/null

    if ! xtb input.xyz --opt normal --gfn 2 --chrg "$chrg" $solv_flag --etemp 500 > opt1.out 2>&1; then
        echo "[ERROR] Stage 1 optimization failed for $name:"
        tail -n 20 opt1.out
        popd > /dev/null
        return 1
    fi

    if ! xtb xtbopt.xyz --opt tight --gfn 2 --chrg "$chrg" $solv_flag > opt2.out 2>&1; then
        cp xtbopt.xyz final.xyz
        cp opt1.out final.out
    else
        cp xtbopt.xyz final.xyz
        cp opt2.out final.out
    fi

    local energy
    energy=$(grep "TOTAL ENERGY" final.out | tail -n 1 | awk '{print $4}')
    echo "[OK] $name -> E = $energy Eh"
    popd > /dev/null
}

run_xtb_opt "dme_gas"   "$WORKDIR/dme_init.xyz"  0 ""
run_xtb_opt "dme_solv"  "$WORKDIR/dme_init.xyz"  0 "--alpb thf"
run_xtb_opt "urea_gas"  "$WORKDIR/urea_init.xyz" 0 ""
run_xtb_opt "urea_solv" "$WORKDIR/urea_init.xyz" 0 "--alpb dmso"

for n in 4 5 6; do
    for lig in dme urea; do
        solv="thf"
        [[ "$lig" == "urea" ]] && solv="dmso"

        run_xtb_opt "mg_${lig}${n}_pre" "$WORKDIR/mg_${lig}${n}_init.xyz" 2 "--alpb $solv"

        mkdir -p "$WORKDIR/mg_${lig}${n}_md"
        cp "$WORKDIR/mg_${lig}${n}_pre/final.xyz" "$WORKDIR/mg_${lig}${n}_md/start.xyz"
        pushd "$WORKDIR/mg_${lig}${n}_md" > /dev/null
        cat << 'MDINP' > md.inp
$md
   temp=350.0
   time=10.0
   dump=1000.0
   step=2.0
   shake=1
$end
$wall
   potential=logfermi
   sphere: 16.0, all
$end
MDINP
        if xtb start.xyz --md --input md.inp --gfn 2 --chrg 2 --alpb "$solv" > md.out 2>&1; then
            natoms=$(head -n 1 start.xyz)
            tail -n $((natoms + 2)) xtb.trj > annealed.xyz
        else
            cp start.xyz annealed.xyz
        fi
        popd > /dev/null

        run_xtb_opt "mg_${lig}${n}_solv" "$WORKDIR/mg_${lig}${n}_md/annealed.xyz" 2 "--alpb $solv"
        run_xtb_opt "mg_${lig}${n}_gas"  "$WORKDIR/mg_${lig}${n}_solv/final.xyz"   2 ""
    done
done

echo "=== Benchmark xTB Runs Complete! ==="
