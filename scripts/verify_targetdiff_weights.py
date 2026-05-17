"""Verify that a freshly-dropped TargetDiff checkpoint is intact and loadable.

Run after copying `pretrained_diffusion.pt` into
`autonomous_drug_discovery/modules/02_generation/targetdiff/pretrained_models/`.

Usage:
    conda run -n targetdiff_env python scripts/verify_targetdiff_weights.py
    # or, from inside the Docker image:
    python scripts/verify_targetdiff_weights.py

Exits non-zero on any failure so it can gate a Makefile / CI step.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
CKPT_PATH = (
    REPO_ROOT
    / "autonomous_drug_discovery"
    / "modules"
    / "02_generation"
    / "targetdiff"
    / "pretrained_models"
    / "pretrained_diffusion.pt"
)
EGNN_PATH = CKPT_PATH.parent / "egnn_pdbbind_v2016.pt"


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while data := f.read(chunk):
            h.update(data)
    return h.hexdigest()


def inspect(path: Path, required: bool) -> int:
    label = path.name
    if not path.exists():
        msg = "MISSING" if required else "missing (optional)"
        print(f"[{msg:>20}] {label}")
        return 1 if required else 0

    size_mb = path.stat().st_size / (1024 * 1024)
    print(f"[{'present':>20}] {label}  {size_mb:.1f} MB  sha256={sha256(path)[:16]}…")

    try:
        import torch
    except ImportError:
        print(
            f"[{'torch unavailable':>20}] cannot deserialize — run inside `targetdiff_env`"
        )
        return 2

    try:
        obj = torch.load(path, map_location="cpu", weights_only=False)
    except Exception as exc:  # noqa: BLE001
        print(f"[{'LOAD FAILED':>20}] {label}: {exc}")
        return 1

    if isinstance(obj, dict):
        keys = list(obj.keys())[:6]
        sd_key = next((k for k in ("model", "state_dict", "model_state_dict") if k in obj), None)
        state = obj[sd_key] if sd_key else obj
        if isinstance(state, dict):
            tensors = [v for v in state.values() if hasattr(v, "numel")]
            total = sum(t.numel() for t in tensors)
            print(
                f"[{'load ok':>20}] top-level keys={keys}  "
                f"{'wrapper=' + sd_key + '  ' if sd_key else ''}"
                f"tensors={len(tensors)}  params≈{total / 1e6:.1f}M"
            )
        else:
            print(f"[{'load ok':>20}] top-level keys={keys}")
    else:
        print(f"[{'load ok':>20}] type={type(obj).__name__}")
    return 0


def main() -> int:
    print(f"Repo root: {REPO_ROOT}")
    print(f"Drop dir:  {CKPT_PATH.parent}")
    print()
    rc = inspect(CKPT_PATH, required=True)
    rc |= inspect(EGNN_PATH, required=False)
    print()
    if rc == 0:
        print("OK — checkpoints look healthy. Next: run a 1-mol CPU smoke sample via")
        print("  conda run -n targetdiff_env python \\")
        print("    autonomous_drug_discovery/modules/02_generation/targetdiff/scripts/sample_for_pocket.py \\")
        print("    autonomous_drug_discovery/modules/02_generation/targetdiff/configs/sampling.yml \\")
        print("    --pdb_path autonomous_drug_discovery/modules/02_generation/targetdiff/examples/1h36_A_rec_1h36_r88_lig_tt_docked_0_pocket10.pdb \\")
        print("    --result_path /tmp/targetdiff_smoke --device cpu --num_samples 1")
    else:
        print("FAILED — inspect output above.")
    return rc


if __name__ == "__main__":
    sys.exit(main())
