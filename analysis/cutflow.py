"""Cutflow booking, extraction, and table formatting."""

import os
import math


def asimov_significance(S, B):
    """Asimov significance: Z_A = sqrt(2 * [(S+B)*ln(1 + S/B) - S])."""
    if B <= 0 or S <= 0:
        return float("nan")
    return math.sqrt(2 * ((S + B) * math.log(1 + S / B) - S))


def book_cutflow_actions(trig_list, cutflow_steps, mc, data_df,
                         mc_acc, data_acc, decay_modes, active_modes,
                         sig_samples, dy_samples, tt_samples, qcd_samples,
                         unified_ptrs):
    """Book cutflow Count/Sum actions into unified_ptrs.

    Returns:
        phase1_data: dict {trig_name: {cutflow_ptrs, data_n_ptr, sum_w, sum_w2, mode_w}}
    """
    _proc_samples = {"DY": dy_samples, "TT": tt_samples,
                     "HHbbtt": sig_samples, "QCD": qcd_samples}
    phase1_data = {}

    for trig_name, trig in trig_list:
        if trig is not None:
            data_cf = data_df.Filter(trig) if data_df is not None else None
            mc_cf = {name: df.Filter(trig) for name, df in mc.items()}
            data_sel = data_acc.Filter(trig) if data_acc is not None else None
            mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}
        else:
            data_cf = data_df
            mc_cf = dict(mc)
            data_sel = data_acc
            mc_sel = dict(mc_acc)

        cutflow_ptrs = []
        for step_name, cut_expr in cutflow_steps:
            if cut_expr is not None:
                if data_cf is not None:
                    data_cf = data_cf.Filter(cut_expr)
                mc_cf = {s: d.Filter(cut_expr) for s, d in mc_cf.items()}
            d_ptr = data_cf.Count() if data_cf is not None else None
            w_ptrs = {s: mc_cf[s].Sum("w") for s in mc_cf}
            cutflow_ptrs.append((step_name, d_ptr, w_ptrs))
            if d_ptr is not None:
                unified_ptrs.append(d_ptr)
            unified_ptrs.extend(w_ptrs.values())

        data_n_ptr = data_sel.Count() if data_sel is not None else None
        sum_w = {s: mc_sel[s].Sum("w") for s in mc_sel}
        sum_w2 = {s: mc_sel[s].Define("w2", "w*w").Sum("w2") for s in mc_sel}

        # Per-mode yield sums
        mode_w = {}
        for mode in active_modes:
            info = decay_modes[mode]
            samples_for_mode = _proc_samples.get(info["process"], [])
            if not samples_for_mode:
                mode_w[mode] = []
                continue
            mode_w[mode] = [mc_sel[s].Filter(f"decayMode=={mode}").Sum("w")
                            for s in samples_for_mode]
            unified_ptrs.extend(mode_w[mode])

        if data_n_ptr is not None:
            unified_ptrs.append(data_n_ptr)
        unified_ptrs.extend(sum_w.values())
        unified_ptrs.extend(sum_w2.values())

        phase1_data[trig_name] = {
            "cutflow_ptrs": cutflow_ptrs,
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
    cutflow_tables = {}

    for trig_name, _trig in trig_list:
        p1 = phase1_data[trig_name]
        print("=" * 63)
        print(f"  {trig_name}")
        print("=" * 63)

        cutflow_rows = []
        for step_name, d_ptr, w_ptrs in p1["cutflow_ptrs"]:
            data_n_cf = d_ptr.GetValue() if d_ptr is not None else -1
            dy_cf = sum(w_ptrs[s].GetValue() for s in dy_samples) * lumi
            tt_cf = sum(w_ptrs[s].GetValue() for s in tt_samples) * lumi
            qcd_cf = sum(w_ptrs[s].GetValue() for s in qcd_samples) * lumi if qcd_samples else 0.0
            sig_cf = sum(w_ptrs[s].GetValue() for s in sig_samples) * lumi
            B_cf = dy_cf + tt_cf + qcd_cf
            s_sqrtb = sig_cf / math.sqrt(B_cf) if B_cf > 0 else float("nan")
            z_a = asimov_significance(sig_cf, B_cf)
            cutflow_rows.append(
                (step_name, data_n_cf, dy_cf, tt_cf, qcd_cf, sig_cf, s_sqrtb, z_a))
        cutflow_tables[trig_name] = cutflow_rows

        print(f"\n{'Cut':<30s} {'Data':>12s} {'DY':>14s} {'TT':>14s} "
              f"{'QCD':>14s} {'Signal':>14s} {'S/√B':>10s} {'Z_A':>10s}")
        print("-" * 122)
        for row in cutflow_rows:
            step, dn, dy, tt, qcd, sig, sb, za = row
            dn_str = f"{dn:>12d}" if dn >= 0 else f"{'(no data)':>12s}"
            print(f"{step:<30s} {dn_str} {dy:>14.4g} {tt:>14.4g} "
                  f"{qcd:>14.4g} {sig:>14.4g} {sb:>10.4g} {za:>10.4g}")

        data_n_ptr = p1["data_n_ptr"]
        sum_w, sum_w2 = p1["sum_w"], p1["sum_w2"]
        mode_w = p1["mode_w"]
        data_n = data_n_ptr.GetValue() if data_n_ptr is not None else -1

        def group_yield(samples):
            y = sum(sum_w[s].GetValue() for s in samples) * lumi
            yerr = math.sqrt(sum(sum_w2[s].GetValue() for s in samples)) * lumi
            return y, yerr

        dy_y, dy_yerr = group_yield(dy_samples)
        tt_y, tt_yerr = group_yield(tt_samples)
        qcd_y, qcd_yerr = group_yield(qcd_samples) if qcd_samples else (0.0, 0.0)
        sig_y, sig_yerr = group_yield(sig_samples)
        B = dy_y + tt_y + qcd_y
        S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")

        print(f"\nData events: {data_n if data_n >= 0 else '(no data)'}")
        print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
        print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
        print(f"QCD yield:   {qcd_y:.6g} +/- {qcd_yerr:.3g}")
        print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")
        print(f"S/sqrt(B):   {S_over_sqrtB:.6g}")
        for mode in active_modes:
            y_mode = sum(d.GetValue() for d in mode_w[mode]) * lumi
            print(f"  mode {mode:>2d} ({decay_modes[mode]['label']:>25s}): {y_mode:.6g}")
        print()

    # Write cutflow markdown
    cutflow_path = os.path.join(plot_dir, "cutflow.md")
    with open(cutflow_path, "w") as f:
        f.write("# Cutflow Tables\n\n")
        f.write(f"Luminosity: {lumi} fb$^{{-1}}$\n\n")
        for trig_name, rows in cutflow_tables.items():
            f.write(f"## {trig_name}\n\n")
            f.write("| Cut | Data | DY | TT | QCD | Signal "
                    "| S/sqrt(B) | Z_A (Asimov) |\n")
            f.write("|-----|-----:|---:|---:|----:|-------:"
                    "|----------:|-------------:|\n")
            for step, dn, dy, tt, qcd, sig, sb, za in rows:
                f.write(f"| {step} | {dn:,d} | {dy:.4g} | {tt:.4g} "
                        f"| {qcd:.4g} | {sig:.4g} | {sb:.4g} | {za:.4g} |\n")
            f.write("\n")
    print(f"Cutflow table saved: {cutflow_path}")

    return cutflow_tables
