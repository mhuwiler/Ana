"""Load analysis configuration from YAML files."""

import os
from dataclasses import dataclass
import yaml

from utils.triggers import DST_JetHT_expr, DST_MU_expr, DST_EL_expr, PARKING_HH_expr

_TRIGGER_EXPR_MAP = {
    "DST_JetHT_expr": DST_JetHT_expr,
    "DST_MU_expr": DST_MU_expr,
    "DST_EL_expr": DST_EL_expr,
    "PARKING_HH_expr": PARKING_HH_expr,
}


@dataclass
class AnalysisConfig:
    """Parsed contents of config/analysis.yaml."""
    decay_modes: dict     # {int: {process, label, color}}
    bkg_modes: list       # list of ints
    sig_modes: list       # list of ints
    cutflow_steps: list   # list of (name, expr_or_None) tuples


def _load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f) or {}


def load_analysis_config(config_dir="config"):
    """Load decay modes and cutflow steps from config/analysis.yaml.

    Returns
    -------
    AnalysisConfig
    """
    cfg = _load_yaml(os.path.join(config_dir, "analysis.yaml"))

    decay_modes_raw = cfg.get("decay_modes", {})
    decay_modes = {int(k): v for k, v in decay_modes_raw.items()}

    bkg_modes = cfg.get("bkg_modes", [1, 2, 3, 4, 5, 10, 11, 12, 30])
    sig_modes = cfg.get("sig_modes", [20, 21, 22])

    cutflow_steps_raw = cfg.get("cutflow_steps", [])
    cutflow_steps = [(s[0], s[1]) for s in cutflow_steps_raw]

    return AnalysisConfig(
        decay_modes=decay_modes,
        bkg_modes=bkg_modes,
        sig_modes=sig_modes,
        cutflow_steps=cutflow_steps,
    )


def load_triggers(config_dir="config"):
    """Load trigger definitions from config/triggers.yaml.

    Returns:
        trig_list: list of (name, expr_string_or_None) tuples
        excl_trig_list: list of (name, expr_string_or_None) tuples
        brilcalc_files: dict {trigger_key: path}
        brilcalc_fallback: dict {trigger_key: fallback_trigger_key}
    """
    cfg = _load_yaml(os.path.join(config_dir, "triggers.yaml"))

    trig_list = []
    for t in cfg.get("triggers", []):
        name = t["name"]
        expr_ref = t.get("expr")
        if expr_ref is None:
            expr = None
        else:
            expr = _TRIGGER_EXPR_MAP.get(expr_ref, expr_ref)
        trig_list.append((name, expr))

    # Build exclusive triggers
    excl_trig_list = []
    if cfg.get("exclusive", True):
        all_exprs = [e for _, e in trig_list if e is not None]
        for name, expr in trig_list:
            if expr is None:
                excl_trig_list.append((name, None))
            else:
                others = " || ".join(e for e in all_exprs if e != expr)
                excl_expr = f"({expr}) && !({others})"
                excl_trig_list.append((f"{name}_excl", excl_expr))

    # Brilcalc file paths
    ana_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    brilcalc_files = {}
    for key, path in cfg.get("brilcalc_files", {}).items():
        brilcalc_files[key] = os.path.join(ana_dir, path)

    brilcalc_fallback = cfg.get("brilcalc_fallback", {})

    return trig_list, excl_trig_list, brilcalc_files, brilcalc_fallback
