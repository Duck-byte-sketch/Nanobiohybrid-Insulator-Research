import json, shutil
from pathlib import Path

v2_path = Path("results/verified_v2_data.json")
fig_dir = Path("poster_figures")
arch_dir = fig_dir / "archive_pending_neb"
arch_dir.mkdir(parents=True, exist_ok=True)

v2 = json.loads(v2_path.read_text(encoding="utf-8"))

for entry in v2.get("vertical_reduction_alpb_thf", {}).values():
    if isinstance(entry, dict):
        sp = entry.get("spin_localization", {})
        if isinstance(sp, dict) and sp.get("mg_spin") == 0 and sp.get("ligand_spin") == 0:
            entry["spin_localization"] = {"status": "not_determined"}
        entry["poster_display"] = False

for k, val in v2.items():
    if ("migration" in k.lower() or "barrier" in k.lower()) and isinstance(val, dict):
        val["poster_status"] = "excluded_pending_neb_validation"

v2_path.write_text(json.dumps(v2, indent=2), encoding="utf-8")

for name in ("fig2_mg_migration_barrier.png", "fig2_mg_migration_barrier.svg"):
    p = fig_dir / name
    if p.exists():
        shutil.move(str(p), str(arch_dir / name))

print("Step 1 complete: JSON sanitized and barrier plots archived.")
