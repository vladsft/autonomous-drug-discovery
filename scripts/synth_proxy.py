#!/usr/bin/env python3
"""Fast synthesizability proxy — the speed tier of the Stage-2.5 gate.

Full AiZynthFinder retrosynthesis costs ~40-50 s/molecule. RAScore is a small
model trained to PREDICT AiZynth solvability from a fingerprint in ~ms/molecule,
so it makes a cheap pre-filter: score everything with the proxy, keep the likely-
makeable, and (optionally) run real AiZynth only on those to obtain actual routes.

Backends, tried in order:
  1. RAScore (reymond-group) — returns P(AiZynth-solvable) in [0,1] from an XGB
     model in ~ms. INSTALL CAVEAT (2026-05-30): the upstream package pins an
     ancient scikit-learn that fails to build on Python 3.11
     (`ModuleNotFoundError: pkg_resources`), so a plain
     `pip install git+https://github.com/reymond-group/RAscore.git` into a
     modern env does NOT work. To activate the fast path, install it in a
     DEDICATED env with the pinned old deps (python 3.8/3.9 + scikit-learn
     ~1.0), mirroring the aizynth_env pattern, and have this module shell out to
     it. Until then the gate runs full AiZynth — correct, just ~30 min instead
     of ~2 for a campaign's survivors.
  2. none — proxy unavailable (the current state). `available()` returns False;
     the gate falls back to running full AiZynth on everything.

We deliberately do NOT ship a hand-rolled heuristic proxy: SA-score-style
heuristics were shown to mis-rate generative scaffolds as "easy" (the very
failure that motivated this gate). Better an honest slow path than a misleading
fast one.
"""
from __future__ import annotations

import functools


@functools.lru_cache(maxsize=1)
def _load_rascore():
    """Return a callable smiles->score, or None if RAScore isn't installed."""
    try:
        from RAscore import RAscore_XGB  # type: ignore
    except Exception:
        return None
    try:
        scorer = RAscore_XGB.RAScorerXGB()
        return scorer.predict
    except Exception as e:  # model load failed
        print(f"[synth_proxy] RAScore present but failed to load ({e}); "
              f"falling back to full AiZynth")
        return None


def available() -> bool:
    """True iff a fast proxy backend is usable."""
    return _load_rascore() is not None


def backend_name() -> str:
    return "rascore" if available() else "none"


def score(smiles_list: list[str]) -> dict[str, float]:
    """Map each SMILES → proxy synthesizability probability in [0,1].

    Returns {} when no proxy is available (caller falls back to full AiZynth).
    Unparseable / erroring SMILES get 0.0 (treated as not-makeable by the proxy,
    but the gate's AiZynth-confirm step is the real arbiter).
    """
    predict = _load_rascore()
    if predict is None:
        return {}
    out: dict[str, float] = {}
    for smi in smiles_list:
        try:
            out[smi] = float(predict(smi))
        except Exception:
            out[smi] = 0.0
    return out
