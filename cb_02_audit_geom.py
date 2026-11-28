import json, math
from pathlib import Path
B = Path("biohybrid_screen_v2")
F = {
  "Mg_DME3": B / "complex_Mg_DME2_DME/xtbopt.xyz",
  "Mg_DME2_MeOEt": B / "complex_Mg_DME2_Methoxyethanol/xtbopt.xyz",
  "Mg_DME2_DMA": B / "complex_Mg_DME2_Amide_DMA/xtbopt.xyz",
  "Mg_DME2_Betaine": B / "complex_Mg_DME2_Zwitterion_Bet/xtbopt.xyz",
  "NEB_pure_start": Path("relaxed_hop_v2/pure_bet_fwd_0/xtbopt.xyz"),
  "NEB_hyb_start": Path("relaxed_hop_v2/hybrid_dma_fwd_0/xtbopt.xyz"),
}
rep = {}
for k, p in F.items():
  if not p.exists(): continue
  ls = p.read_text().strip().splitlines()
  n = int(ls[0].strip())
  ats = [(l.split()[0], [float(x) for x in l.split()[1:4]]) for l in ls[2:2+n]]
  mgs = [c for s, c in ats if s == "Mg"]
  bs = sorted([round(math.dist(mgs[0], c), 3) for s, c in ats if s == "O" and math.dist(mgs[0], c) < 2.6]) if mgs else []
  mc = [i for i, (s, c) in enumerate(ats) if s == "C" and sum(1 for s2, c2 in ats if s2 == "H" and math.dist(c, c2) < 1.25) == 3]
  rep[k] = {"n_atoms": n, "Mg_O_A": bs, "CN": len(bs), "methyl_C": mc}
Path("results").mkdir(exist_ok=True)
Path("results/v3_geom_audit.json").write_text(json.dumps(rep, indent=2))
print("SAVED_AUDIT:", list(rep.keys()))
