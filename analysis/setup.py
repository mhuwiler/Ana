"""ROOT initialization and C++ macro loading."""

import os
import resource

from analysis.constants import STACK_SIZE_MB


def init_root(nthreads=4, stack_mb=None):
    """Initialize ROOT with ImplicitMT and increased stack size.

    Parameters
    ----------
    nthreads : int
        Number of threads for ImplicitMT.
    stack_mb : int or None
        Stack size in MB. Defaults to STACK_SIZE_MB from constants.
    """
    if stack_mb is None:
        stack_mb = STACK_SIZE_MB
    stack_bytes = stack_mb * 1024 * 1024
    resource.setrlimit(resource.RLIMIT_STACK, (stack_bytes, resource.RLIM_INFINITY))

    import ROOT
    ROOT.gErrorIgnoreLevel = ROOT.kInfo
    ROOT.ROOT.EnableImplicitMT(nthreads)


def load_macros(ana_dir=None):
    """Load C++ analysis macros. Uses ACLiC if .so exists, else interpreted.

    Parameters
    ----------
    ana_dir : str or None
        Root directory of the analysis repo. Defaults to parent of this file's
        directory.
    """
    import ROOT

    if ana_dir is None:
        ana_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    os.chdir(ana_dir)
    import sys
    sys.path.insert(0, ana_dir)
    ROOT.gInterpreter.AddIncludePath(ana_dir)
    ROOT.gInterpreter.AddIncludePath(os.path.join(ana_dir, "elements"))

    for macro in ["elements/GenMatching.C", "elements/RecoObjects.C"]:
        so = macro.replace(".C", "_C.so")
        if os.path.exists(so) and os.path.getmtime(so) >= os.path.getmtime(macro):
            ROOT.gROOT.LoadMacro(f"{macro}+")
        else:
            ROOT.gROOT.LoadMacro(macro)
    print("CWD =", os.getcwd())
