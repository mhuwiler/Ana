"""Cutflow booking, extraction, and table formatting."""

import os
import math


def asimov_significance(S, B):
    """Asimov significance: Z_A = sqrt(2 * [(S+B)*ln(1 + S/B) - S])."""
    if B <= 0 or S <= 0:
        return float("nan")
    val = 2 * ((S + B) * math.log(1 + S / B) - S)
    if val < 0:
        return float("nan")
    return math.sqrt(val)


def book_cutflow_actions(trig_list, cuts_by_trig, mc, data_df,
                         mc_acc, data_acc, decay_modes, active_modes,
                         sig_samples, dy_samples, tt_samples, qcd_samples,
                         unified_ptrs):
    """Book cutflow Count/Sum actions into unified_ptrs.

    Parameters
    ----------
    cuts_by_trig : dict {trig_name: [(cut_name, expr), ...]}
        Per-trigger cut lists. Cuts are applied progressively for cutflow.

    Returns:
        phase1_data: dict {trig_name: {cutflow_ptrs, data_n_ptr, sum_w, sum_w2, mode_w}}
    """
    _proc_samples = {"DY": dy_samples, "TT": tt_samples,
                     "HHbbtt": sig_samples, "QCD": qcd_samples}
    phase1_data = {}

    for trig_name, trig in trig_list:
        trig_cuts = cuts_by_trig.get(trig_name, [])
        if not trig_cuts:
            continue

        # Build cutflow steps: Trigger first, then each cut
        # Use expression as step name if no name provided
        named_cuts = [(f"{name}: {expr}" if name else expr, expr) for name, expr in trig_cuts]
        cutflow_steps = [("Trigger", None)] + named_cuts

        # Start from pre-cut mc, apply trigger
        if trig is not None:
            data_cf = data_df.Filter(trig) if data_df is not None else None
            mc_cf = {name: df.Filter(trig) for name, df in mc.items()}
        else:
            data_cf = data_df
            mc_cf = dict(mc)

        # Signal modes for per-channel S/√B
        sig_modes = [m for m in active_modes if decay_modes[m]["process"] == "HHbbtt"]

        cutflow_ptrs = []
        for step_name, cut_expr in cutflow_steps:
            if cut_expr is not None:
                if data_cf is not None:
                    data_cf = data_cf.Filter(cut_expr)
                mc_cf = {s: d.Filter(cut_expr) for s, d in mc_cf.items()}
            d_ptr = data_cf.Count() if data_cf is not None else None
            w_ptrs = {s: mc_cf[s].Sum("w") for s in mc_cf}
            # Per-signal-mode yields at this cutflow step
            sig_mode_ptrs = {}
            for mode in sig_modes:
                sig_mode_ptrs[mode] = [mc_cf[s].Filter(f"decayMode=={mode}").Sum("w")
                                       for s in sig_samples]
                unified_ptrs.extend(sig_mode_ptrs[mode])
            cutflow_ptrs.append((step_name, d_ptr, w_ptrs, sig_mode_ptrs))
            if d_ptr is not None:
                unified_ptrs.append(d_ptr)
            unified_ptrs.extend(w_ptrs.values())

        # Final yields (after all cuts) — use mc_cf which has all cuts applied
        data_n_ptr = data_cf.Count() if data_cf is not None else None
        sum_w = {s: mc_cf[s].Sum("w") for s in mc_cf}
        sum_w2 = {s: mc_cf[s].Define("w2", "w*w").Sum("w2") for s in mc_cf}

        # Per-mode yield sums (after all cuts)
        mode_w = {}
        for mode in active_modes:
            info = decay_modes[mode]
            samples_for_mode = _proc_samples.get(info["process"], [])
            if not samples_for_mode:
                mode_w[mode] = []
                continue
            mode_w[mode] = [mc_cf[s].Filter(f"decayMode=={mode}").Sum("w")
                            for s in samples_for_mode]
            unified_ptrs.extend(mode_w[mode])

        if data_n_ptr is not None:
            unified_ptrs.append(data_n_ptr)
        unified_ptrs.extend(sum_w.values())
        unified_ptrs.extend(sum_w2.values())

        phase1_data[trig_name] = {
            "cutflow_ptrs": cutflow_ptrs,
            "cutflow_cards": [c for _, c_list in [(trig_name, trig_cuts)]
                              for c in ["cuts"]],
            "data_n_ptr": data_n_ptr,
            "sum_w": sum_w, "sum_w2": sum_w2,
            "mode_w": mode_w,
        }

    return phase1_data


def extract_and_print_cutflow(trig_list, phase1_data, decay_modes, active_modes,
                              dy_samples, tt_samples, qcd_samples, sig_samples,
                              lumi, plot_dir):
    """Extract cutflow results and print + save markdown table.

    Returns:
        cutflow_tables: dict {trig_name: rows}
    """
    sig_modes = [m for m in active_modes if decay_modes[m]["process"] == "HHbbtt"]
    sig_labels = {m: decay_modes[m]["label"] for m in sig_modes}
    # Short labels for column headers
    sig_short = {20: "τhτh", 21: "τμτh", 22: "τeτh"}

    cutflow_tables = {}

    for trig_name, _trig in trig_list:
        if trig_name not in phase1_data:
            continue
        p1 = phase1_data[trig_name]
        print("=" * 80)
        print(f"  {trig_name}")
        print("=" * 80)

        cutflow_rows = []
        for step_name, d_ptr, w_ptrs, sig_mode_ptrs in p1["cutflow_ptrs"]:
            data_n_cf = d_ptr.GetValue() if d_ptr is not None else -1
            dy_cf = sum(w_ptrs[s].GetValue() for s in dy_samples) * lumi
            tt_cf = sum(w_ptrs[s].GetValue() for s in tt_samples) * lumi
            qcd_cf = sum(w_ptrs[s].GetValue() for s in qcd_samples) * lumi if qcd_samples else 0.0
            B_cf = max(dy_cf + tt_cf + qcd_cf, 0.0)

            # Per-channel signal yields
            sig_by_mode = {}
            for mode in sig_modes:
                sig_by_mode[mode] = sum(p.GetValue() for p in sig_mode_ptrs[mode]) * lumi

            sig_total = sum(sig_by_mode.values())
            s_sqrtb_total = sig_total / math.sqrt(B_cf) if B_cf > 0 else float("nan")
            s_sqrtb_by_mode = {m: sig_by_mode[m] / math.sqrt(B_cf) if B_cf > 0
                               else float("nan") for m in sig_modes}

            cutflow_rows.append({
                "step": step_name, "data": data_n_cf,
                "dy": dy_cf, "tt": tt_cf, "qcd": qcd_cf, "B": B_cf,
                "sig_by_mode": sig_by_mode, "sig_total": sig_total,
                "s_sqrtb_total": s_sqrtb_total, "s_sqrtb_by_mode": s_sqrtb_by_mode,
            })
        cutflow_tables[trig_name] = cutflow_rows

        # Print table header
        mode_headers = "".join(f" {'S(' + sig_short.get(m, str(m)) + ')':>10s}" for m in sig_modes)
        sqrtb_headers = "".join(f" {'S/√B(' + sig_short.get(m, str(m)) + ')':>12s}" for m in sig_modes)
        print(f"\n{'Cut':<22s} {'Data':>10s} {'DY':>12s} {'TT':>12s} "
              f"{'QCD':>12s} {'B':>12s}{mode_headers}{sqrtb_headers}")
        print("-" * (22 + 10 + 12*4 + (10 + 12) * len(sig_modes) + len(sig_modes) * 2 + 5))

        for row in cutflow_rows:
            dn_str = f"{row['data']:>10d}" if row['data'] >= 0 else f"{'(no data)':>10s}"
            mode_vals = "".join(f" {row['sig_by_mode'][m]:>10.4g}" for m in sig_modes)
            sqrtb_vals = "".join(f" {row['s_sqrtb_by_mode'][m]:>12.4g}" for m in sig_modes)
            print(f"{row['step']:<22s} {dn_str} {row['dy']:>12.4g} {row['tt']:>12.4g} "
                  f"{row['qcd']:>12.4g} {row['B']:>12.4g}{mode_vals}{sqrtb_vals}")

        # Final summary
        data_n_ptr = p1["data_n_ptr"]
        sum_w, sum_w2 = p1["sum_w"], p1["sum_w2"]
        mode_w = p1["mode_w"]
        data_n = data_n_ptr.GetValue() if data_n_ptr is not None else -1

        def group_yield(samples):
            y = sum(sum_w[s].GetValue() for s in samples) * lumi
            yerr = math.sqrt(max(sum(sum_w2[s].GetValue() for s in samples), 0.0)) * lumi
            return y, yerr

        dy_y, dy_yerr = group_yield(dy_samples)
        tt_y, tt_yerr = group_yield(tt_samples)
        qcd_y, qcd_yerr = group_yield(qcd_samples) if qcd_samples else (0.0, 0.0)
        sig_y, sig_yerr = group_yield(sig_samples)

        print(f"\nFinal yields (after all cuts):")
        print(f"  DY:     {dy_y:>12.4g} +/- {dy_yerr:.3g}")
        print(f"  TT:     {tt_y:>12.4g} +/- {tt_yerr:.3g}")
        print(f"  QCD:    {qcd_y:>12.4g} +/- {qcd_yerr:.3g}")
        print(f"  Signal: {sig_y:>12.4g} +/- {sig_yerr:.3g}")
        for mode in sig_modes:
            y_mode = sum(d.GetValue() for d in mode_w[mode]) * lumi
            print(f"    {sig_labels[mode]:>30s}: {y_mode:.4g}")
        print()

    # Write cutflow markdown
    cutflow_path = os.path.join(plot_dir, "cutflow.md")
    with open(cutflow_path, "w") as f:
        f.write("# Cutflow Tables\n\n")
        f.write(f"Luminosity: {lumi} fb$^{{-1}}$\n\n")
        for trig_name, rows in cutflow_tables.items():
            f.write(f"## {trig_name}\n\n")
            # Header
            cols = "| Cut | DY | TT | QCD | B |"
            for m in sig_modes:
                cols += f" S({sig_short.get(m, str(m))}) |"
            for m in sig_modes:
                cols += f" S/√B({sig_short.get(m, str(m))}) |"
            f.write(cols + "\n")
            f.write("|" + "|".join(["-----"] + ["---:"] * (4 + 2 * len(sig_modes))) + "|\n")
            for row in rows:
                line = f"| {row['step']} | {row['dy']:.4g} | {row['tt']:.4g} "
                line += f"| {row['qcd']:.4g} | {row['B']:.4g} |"
                for m in sig_modes:
                    line += f" {row['sig_by_mode'][m]:.4g} |"
                for m in sig_modes:
                    line += f" {row['s_sqrtb_by_mode'][m]:.4g} |"
                f.write(line + "\n")
            f.write("\n")
    print(f"Cutflow table saved: {cutflow_path}")

    return cutflow_tables
