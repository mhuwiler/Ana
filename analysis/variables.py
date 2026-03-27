"""PlotVar dataclass and helper functions for plot variable management.

This module provides the PlotVar and PlotVar2D types and expansion helpers.
Variable *instances* are defined in user scripts, not here.
"""

import fnmatch
from dataclasses import dataclass

from analysis.constants import AK8_PT_MIN


@dataclass(frozen=True)
class PlotVar:
    """A single 1D plot variable definition."""
    name: str
    xlabel: str
    nbins: int
    xmin: float
    xmax: float
    sentinel: bool = False  # True for AK8/gen-matched vars with -1 sentinel


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


def expand_ak8_vars(count=3, labels=None, pt_maxes=None, vars_per_jet=None):
    """Generate per-jet AK8 PlotVars from a template.

    Parameters
    ----------
    count : int
        Number of leading AK8 jets.
    labels : list[str]
        Human-readable label per jet (e.g. ["Leading", "Sub-leading", "Third"]).
    pt_maxes : list[float]
        pT axis upper bound per jet.
    vars_per_jet : list[tuple]
        Each tuple: (suffix, xlabel_template, nbins, xmin, xmax_or_marker).
        Use "{pt_max}" as xmax to substitute the per-jet pt_max.

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
            full_name = f"ak8_{suffix}{i}"
            full_xlabel = f"{lbl} AK8 {xlabel_tmpl}"
            xmax = pt_max if str(xmax_tmpl) == "{pt_max}" else xmax_tmpl
            xmin_val = AK8_PT_MIN if suffix == "pt" else xmin
            result.append(PlotVar(full_name, full_xlabel, nbins, xmin_val, xmax,
                                  sentinel=True))
    return result


def expand_ak8_candidate_vars(candidates, candidate_vars):
    """Generate AK8 candidate PlotVars (e.g. Hbb, Htautau candidates).

    Parameters
    ----------
    candidates : list[tuple[str, str]]
        Each tuple: (candidate_name, candidate_label).
    candidate_vars : list[tuple]
        Each tuple: (suffix, xlabel, nbins, xmin, xmax).

    Returns
    -------
    list[PlotVar]
    """
    result = []
    for cand_name, cand_label in candidates:
        for var_def in candidate_vars:
            suffix, xlabel, nbins, xmin, xmax = var_def
            full_name = f"ak8_{cand_name}_{suffix}"
            full_xlabel = f"{cand_label} cand AK8 {xlabel}"
            result.append(PlotVar(full_name, full_xlabel, nbins, xmin, xmax,
                                  sentinel=True))
    return result


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
