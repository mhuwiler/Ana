"""PlotVar dataclass and helper functions for plot variable management.

This module provides the PlotVar and PlotVar2D types and expansion helpers.
Variable *instances* are defined in user scripts, not here.
"""

import fnmatch
from dataclasses import dataclass

from analysis.constants import AK8_PT_MIN


def _auto_sentinel(name):
    """Infer whether a variable uses -1 sentinel from its name prefix."""
    return name.startswith("ak8_") or name.startswith("ak4_gen") or name.startswith("gen_")


@dataclass(frozen=True)
class PlotVar:
    """A single 1D plot variable definition."""
    name: str
    xlabel: str
    nbins: int
    xmin: float
    xmax: float
    sentinel: bool = None  # None = auto-detect from name prefix


@dataclass(frozen=True)
class PlotVar2D:
    """A single 2D plot variable definition."""
    name: str
    xvar: str
    yvar: str
    xlabel: str
    ylabel: str
    nx: int
    x0: float
    x1: float
    ny: int
    y0: float
    y1: float


def expand_jet_vars(prefix, count=3, labels=None, pt_maxes=None,
                    vars_per_jet=None, pt_min=0, sentinel=True):
    """Generate per-jet PlotVars from a template.

    Parameters
    ----------
    prefix : str
        Variable name prefix (e.g. "ak8", "ak4").
    count : int
        Number of leading jets.
    labels : list[str]
        Human-readable label per jet (e.g. ["Leading", "Sub-leading", "Third"]).
    pt_maxes : list[float]
        pT axis upper bound per jet.
    vars_per_jet : list[tuple]
        Each tuple: (suffix, xlabel_template, nbins, xmin, xmax_or_marker).
        Use "{pt_max}" as xmax to substitute the per-jet pt_max.
    pt_min : float
        Override xmin for "pt" suffix variables.
    sentinel : bool
        Whether these variables use -1 sentinel filtering.

    Returns
    -------
    list[PlotVar]
    """
    if labels is None:
        labels = [f"Jet {i}" for i in range(count)]
    if pt_maxes is None:
        pt_maxes = [1000] * count
    if vars_per_jet is None:
        return []

    result = []
    for i in range(count):
        lbl = labels[i] if i < len(labels) else f"Jet {i}"
        pt_max = pt_maxes[i] if i < len(pt_maxes) else 1000
        for var_def in vars_per_jet:
            suffix, xlabel_tmpl, nbins, xmin, xmax_tmpl = var_def
            full_name = f"{prefix}_{suffix}{i}"
            full_xlabel = f"{lbl} {prefix.upper()} {xlabel_tmpl}"
            xmax = pt_max if str(xmax_tmpl) == "{pt_max}" else xmax_tmpl
            xmin_val = pt_min if suffix == "pt" and pt_min > 0 else xmin
            result.append(PlotVar(full_name, full_xlabel, nbins, xmin_val, xmax,
                                  sentinel=sentinel))
    return result


def expand_candidate_vars(prefix, candidates, candidate_vars, sentinel=True):
    """Generate candidate PlotVars (e.g. Hbb, Htautau candidates).

    Parameters
    ----------
    prefix : str
        Variable name prefix (e.g. "ak8").
    candidates : list[tuple[str, str]]
        Each tuple: (candidate_name, candidate_label).
    candidate_vars : list[tuple]
        Each tuple: (suffix, xlabel, nbins, xmin, xmax).
    sentinel : bool
        Whether these variables use -1 sentinel filtering.

    Returns
    -------
    list[PlotVar]
    """
    result = []
    for cand_name, cand_label in candidates:
        for var_def in candidate_vars:
            suffix, xlabel, nbins, xmin, xmax = var_def
            full_name = f"{prefix}_{cand_name}_{suffix}"
            full_xlabel = f"{cand_label} cand {prefix.upper()} {xlabel}"
            result.append(PlotVar(full_name, full_xlabel, nbins, xmin, xmax,
                                  sentinel=sentinel))
    return result


def expand_ak8_vars(count=3, labels=None, pt_maxes=None, vars_per_jet=None):
    """Generate per-jet AK8 PlotVars. Convenience wrapper for expand_jet_vars."""
    return expand_jet_vars("ak8", count=count, labels=labels, pt_maxes=pt_maxes,
                           vars_per_jet=vars_per_jet, pt_min=AK8_PT_MIN, sentinel=True)


def expand_ak8_candidate_vars(candidates, candidate_vars):
    """Generate AK8 candidate PlotVars. Convenience wrapper for expand_candidate_vars."""
    return expand_candidate_vars("ak8", candidates, candidate_vars, sentinel=True)


def var2d_attrs(v):
    """Return (name, xvar, yvar, xlabel, ylabel, nx, x0, x1, ny, y0, y1) from a PlotVar2D or tuple."""
    if isinstance(v, PlotVar2D):
        return v.name, v.xvar, v.yvar, v.xlabel, v.ylabel, v.nx, v.x0, v.x1, v.ny, v.y0, v.y1
    return tuple(v)


def var_attrs(v):
    """Return (name, xlabel, nbins, xmin, xmax) from a PlotVar or 5-tuple.

    Allows histogramming and plotting code to accept either a PlotVar instance
    or a plain tuple during the transition period.
    """
    if isinstance(v, PlotVar):
        return v.name, v.xlabel, v.nbins, v.xmin, v.xmax
    return v[0], v[1], v[2], v[3], v[4]


def is_sentinel(v):
    """Return True if this variable uses -1 sentinel filtering (AK8 / gen-matched).

    Accepts PlotVar (uses .sentinel field, auto-detects if None) or tuple (infers from name).
    """
    if isinstance(v, PlotVar):
        return _auto_sentinel(v.name) if v.sentinel is None else v.sentinel
    name = v[0]
    return _auto_sentinel(name)


def filter_vars(vars_list, patterns):
    """Apply fnmatch filter for --plot-vars CLI flag.

    Returns all vars if patterns is None.
    """
    if patterns is None:
        return vars_list
    filtered = [v for v in vars_list
                if any(fnmatch.fnmatch(v.name, p) for p in patterns)]
    if not filtered:
        all_names = {v.name for v in vars_list}
        unknown = {p for p in patterns if not any(fnmatch.fnmatch(n, p) for n in all_names)}
        if unknown:
            print(f"WARNING: no variables matched patterns: {unknown}")
    return filtered
