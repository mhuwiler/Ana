#!/usr/bin/env python3
"""Collect analysis plots into a PDF.

Auto-discovers from output/plots/{script}/ directory. One page per
(variable, plot_type) with triggers arranged in a grid.

Usage:
    # Stacked + sig plots from coi_study
    python make_pdf.py --dir output/plots/coi_study --types stacked sig

    # Filter variables
    python make_pdf.py --dir output/plots/coi_study --vars "b_coi*" "tau_coi*"

    # Specific triggers
    python make_pdf.py --dir output/plots/coi_study --triggers NoTrigger DST_JetHT

    # All from allMC_wocuts
    python make_pdf.py --dir output/plots/allMC_wocuts --types stacked shape sig

    # Custom output
    python make_pdf.py --dir output/plots/coi_study -o my_report.pdf
"""

import argparse
import os

import importlib.util, sys
_spec = importlib.util.spec_from_file_location("pdf", "utils/pdf.py")
_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_mod)
PDFReport = _mod.PDFReport

parser = argparse.ArgumentParser(description="Collect analysis plots into a PDF")
parser.add_argument("--dir", required=True,
                    help="Plot directory (e.g. output/plots/coi_study)")
parser.add_argument("-o", "--output", default=None,
                    help="Output PDF path (default: {dir}/report.pdf)")
parser.add_argument("--types", nargs="+", default=["stacked"],
                    help="Plot types as separate pages (default: stacked)")
parser.add_argument("--vars", nargs="+", default=None, metavar="PATTERN",
                    help="Variable patterns (fnmatch, e.g. 'b_coi*')")
parser.add_argument("--triggers", nargs="+", default=None,
                    help="Trigger names (default: all found)")
parser.add_argument("--cols", type=int, default=2,
                    help="Grid columns per page (default: 2)")
args = parser.parse_args()

pdf = PDFReport(args.dir)

var_patterns = args.vars or ["*"]
variables = pdf.variables(var_patterns)

if not variables:
    print(f"No variables found matching {var_patterns} in {args.dir}")
    raise SystemExit(1)

print(f"[make_pdf] {len(variables)} variables, "
      f"{len(args.triggers or pdf.triggers())} triggers, "
      f"{len(args.types)} plot types")

for var in variables:
    for ptype in args.types:
        pdf.page_grid(var, ptype, triggers=args.triggers, cols=args.cols)

script_name = os.path.basename(os.path.normpath(args.dir))
output = args.output or os.path.join(args.dir, f"{script_name}.pdf")
pdf.save(output)
