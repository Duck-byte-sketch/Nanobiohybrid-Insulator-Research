
Nanobiohybrid Insulator Design for Magnesium-Ion Batteries

Computational R&E Project (GFN2-xTB / ALPB-THF & Direct-SCF B3LYP-D3BJ/def2-SVP & def2-TZVP/IEF-PCM)

Overview

Multivalent $\text{Mg}^{2+}$ batteries suffer from sluggish desolvation of tightly chelated ether electrolytes (such as $\text{[Mg(DME)}_3\text{]}^{2+}$) and parasitic passivation at the Mg-metal anode. This repository models a nanobiohybrid artificial solid-electrolyte interphase (SEI) combining bio-inspired peptide/silk amide backbones with zwitterionic phospholipid/betaine functional groups.

Key Results Summary


Phase 1: Solvation Shell & Method Calibration (AUDIT_RECORD.txt, pyscf_validate.py)

• Coordination Numbers: $\text{[Mg(DME)}_n\text{]}^{2+}$ ($n=4,5,6$) maintains an octahedral inner shell ($\text{CN}=6$, $\text{Mg–O} = 1.93\text{–}2.03\text{ \AA}$) with outer-shell DMEs at $3.77\text{–}4.19\text{ \AA}$. Monodentate $\text{[Mg(Urea)}_n\text{]}^{2+}$ adopts $\text{CN}=4,5,6$.
• Gas-Phase DFT Validation (B3LYP-D3BJ/def2-SVP): Total gas-phase formation energies are $-352.57\text{ kcal/mol}$ for $\text{[Mg(DME)}_4\text{]}^{2+}$ and $-196.70\text{ kcal/mol}$ for $\text{[Mg(Urea)}_4\text{]}^{2+}$. GFN2-xTB reproduces the qualitative ligand-exchange preference ($+84.33$ vs. $+155.88\text{ kcal/mol}$).

Phase 2 & 3 (V2 Verified): Solvated DME-Displacement Screening (BIOHYBRID_SCREEN.txt, V2_COMPLETE_BENCHMARK.txt, results/verified_v2_data.json)

Reaction in implicit THF: $[\text{Mg(DME)}_3]^{2+} + \text{BioMotif} \longrightarrow [\text{Mg(DME)}_2(\text{BioMotif})]^{2+} + \text{DME}$

Functional Motif	Bio-Polymer Analogue	CN	$\Delta E_{\text{swap}}$ (GFN2-xTB / ALPB-THF)	$\Delta E_{\text{swap}}$ (B3LYP-D3BJ/def2-TZVP / IEF-PCM-THF)	Role
DME (Reference)	Baseline Ether Solvent	6	$0.00\text{ kcal/mol}$	$0.00\text{ kcal/mol}$	Symmetric Tris-Chelate
Alkoxy-OH (MeOEtOH)	Cellulose / Chitosan	5	$+23.41\text{ kcal/mol}$	$+19.02\text{ kcal/mol}$	Endothermic; too weak to strip DME
Amide (DMA)	Silk Fibroin / Peptide	5	$+10.77\text{ kcal/mol}$	$+9.66\text{ kcal/mol}$	Endothermic alone; acts as bulk hopping bridge
Zwitterion (Betaine)	Phospholipid / Betaine	6	$-12.35\text{ kcal/mol}$	$-9.57\text{ kcal/mol}$	Sole exothermic interfacial DME stripper

• Solvated DFT Consistency: Upgrading from gas-phase B3LYP-D3BJ/def2-SVP ($-36.92\text{ kcal/mol}$) to implicit-solvent B3LYP-D3BJ/def2-TZVP/IEF-PCM(THF) ($-9.57\text{ kcal/mol}$) removes gas-phase electrostatic overbinding and agrees within ~3–4 kcal/mol of GFN2-xTB/ALPB(THF) ($-12.35\text{ kcal/mol}$) across all tested motifs.

Phase 4 & V2 Diagnostics: Constrained $\text{Mg}^{2+}$ Hopping & Pending UKS Reduction Benchmarks

• Constrained $\text{Mg}^{2+}$ Migration Scans (phase4b_relaxed_scan.py, results/verified_v2_data.json):
	• Single-branch constrained scans with frozen $\text{Mg}^{2+}$ coordinates (atom 1 in $fix) exhibit strong path dependence between forward (-1.80 → +1.80 Å) and reverse (+1.80 → -1.80 Å) branches due to host conformer hysteresis (e.g., forward barriers of $46.65\text{ kcal/mol}$ for pure Betaine vs. $25.33\text{ kcal/mol}$ with an intermediate DMA pocket, whereas the reverse hybrid scan drops into a central pocket).
	• Status: Single-branch constrained scans are path-dependent and barriers are not yet converged; phase4b_relaxed_scan.py implements a 4-seed multi-start lowest-energy envelope using absolute Hartree energies and xtb --path endpoint relaxation to resolve whether the amide pocket acts as a kinetic bridge or a thermodynamic trap.
• Reduction Stability (reduction_alpb_thf/ — Pending UKS/ma-def2-TZVP/IEF-PCM Validation):
	• Preliminary GFN2-xTB/ALPB(THF) vertical electron affinities (VEAs) overestimate the reducibility of neutral amides ($\text{DMA VEA} = +3.56\text{ eV}$) and have unpopulated spin-density entries (0.000). Accurate cathodic stability limits and unpaired-electron localization require unrestricted DFT (UKS B3LYP-D3BJ/ma-def2-TZVP/IEF-PCM(THF)) calibrated against known electrolyte reference molecules.

Repository Structure

• audit.py, build_and_run.py, fix_dme_and_parse.py, AUDIT_RECORD.txt: Phase 1 GFN2-xTB coordination and stepwise desolvation benchmarks.
• pyscf_validate.py, pyscf_verify_betaine.py: Memory-safe Direct-SCF B3LYP-D3BJ validation scripts (def2-SVP and def2-TZVP/IEF-PCM).
• biohybrid_screen_v2/, BIOHYBRID_SCREEN.txt, V2_COMPLETE_BENCHMARK.txt: V2 solvated ligand-displacement optimizations and TZVP/IEF-PCM single-point benchmarks.
• phase4b_relaxed_scan.py: 4-seed multi-start lowest-energy envelope scanner and xtb --path reaction-path runner storing absolute Hartree energies.
• results/verified_v2_data.json: Consolidated JSON archive of verified V2 swap energies and diagnostic scan profiles.
• make_poster_figures.py, poster_figures/: Publication-grade 300-DPI PNG and SVG figure generator and output graphics.
