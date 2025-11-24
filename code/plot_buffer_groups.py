#!/usr/bin/env python3
"""
plot_buffer_groups.py
=====================

Plot buffer base/acid species from UKF outputs, with optional HH overlays.

Inputs:
  reports/softsensor_run/
    - ivt_ukf_results_egfp_HEPES.csv
    - ivt_ukf_results_egfp_TRIS.csv
    - ivt_ukf_results_csp_HEPES.csv
    - ivt_ukf_results_csp_TRIS.csv
Optional HH workbook:
  data/H_H_NTP.xlsx  (sheets named: egfphepes, egfptris, csphepes, csptris)

Outputs:
  figures/buffers/buffers_ALL.pdf

Author: Mahdi Ahmed (2025)
License: MIT
"""
from __future__ import annotations

from pathlib import Path
import argparse
import re

import numpy as np
import pandas as pd
import matplotlib
if not matplotlib.get_backend().lower().startswith("qt"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

AX_LABEL_SIZE    = 10
TICK_LABEL_SIZE  = 10
PANEL_LABEL_SIZE = 10
LINEWIDTH        = 1.8

FIGSIZE       = (7.3, 7.3)
TOP_FIG       = 0.70   # ↓ smaller number = MORE headroom for legend
BOT_FIG       = 0.12
LEFT_FIG      = 0.11
RIGHT_FIG     = 0.995
WSPACE        = 0.38
HSPACE        = 0.50

LEGEND_Y      = 0.85
LEGEND_FS     = 10

PANEL_TITLE_Y = 1.015
TITLE_X_OFF   = 0.10

XLABEL_PAD    = 2
YLABEL_PAD    = 2
XTICK_PAD     = 2
YTICK_PAD     = 2

RUN_LABELS = {
    "egfphepes": "eGFP + HEPES",
    "egfptris":  "eGFP + TRIS",
    "csphepes":  "CSP + HEPES",
    "csptris":   "CSP + TRIS",
}

CSV_DEFAULTS = {
    "eGFP + HEPES": "ivt_ukf_results_egfp_HEPES.csv",
    "eGFP + TRIS":  "ivt_ukf_results_egfp_TRIS.csv",
    "CSP + HEPES":  "ivt_ukf_results_csp_HEPES.csv",
    "CSP + TRIS":   "ivt_ukf_results_csp_TRIS.csv",
}

def pick_column(df, include_keys=(), exclude_keys=()):
    """Looser column matcher: tokenizes and checks includes/excludes."""
    cols = list(df.columns)
    include_keys = list(include_keys)
    exclude_keys = list(exclude_keys)
    for col in cols:
        lc = col.lower()
        tokens = re.split(r'[^a-z0-9]+', lc)
        tokens_set = set(t for t in tokens if t)
        include_ok = True
        for k in include_keys:
            subkeys = [sub for sub in re.split(r'[^a-z0-9]+', k.lower()) if sub]
            if not all(sub in tokens_set for sub in subkeys):
                include_ok = False
                break
        if not include_ok:
            continue
        exclude_flag = False
        for ex in exclude_keys:
            sub_exs = [sub for sub in re.split(r'[^a-z0-9]+', ex.lower()) if sub]
            if any(sub in tokens_set for sub in sub_exs):
                exclude_flag = True
                break
        if exclude_flag:
            continue
        return col
    return None

def _normalize_name(s: str) -> str:
    return re.sub(r'[^a-z0-9]+', '', s.lower())

def _apply_axis_style(ax, y_unit: str):
    ax.set_xlabel('Time (minutes)', fontsize=AX_LABEL_SIZE, labelpad=XLABEL_PAD)
    ax.set_ylabel(f'Concentration ({y_unit})', fontsize=AX_LABEL_SIZE, labelpad=YLABEL_PAD)
    ax.set_xlim(0, 120)
    ax.set_xticks(np.arange(0, 121, 20))
    ax.tick_params(axis='x', pad=XTICK_PAD, labelsize=TICK_LABEL_SIZE)
    ax.tick_params(axis='y', pad=YTICK_PAD, labelsize=TICK_LABEL_SIZE)
    ax.ticklabel_format(style='plain', axis='y', useOffset=False)
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(True)

def _short_label(label: str) -> str:
    nm = " ".join(label.replace('+','').split())
    nm = re.sub(r'\bHEPES\b','',nm,flags=re.IGNORECASE)
    nm = re.sub(r'\bTRIS\b','',nm,flags=re.IGNORECASE)
    return " ".join(nm.split())

def _sheet_map_from_excel(xls: pd.ExcelFile) -> dict[str, str]:
    return {_normalize_name(n): n for n in xls.sheet_names}

def plot_buffer_groups(
    dataset_paths: dict[str, str | Path],
    buffer_groups: list[dict] | None = None,
    smooth: bool = False,
    window: int = 5,
    cols: int = 2,
    figsize: tuple = FIGSIZE,
    split_by_tris_hepes: bool = False,
    hhpred_excel: str | Path | None = None,
    out_pdf: str | Path = "figures/buffers/buffers_ALL.pdf",
):
    """
    dataset_paths: mapping like {"eGFP + HEPES": "reports/...csv", ...}
    buffer_groups: list of dicts defining base/acid field names and units
    hhpred_excel: optional Excel to overlay HH predictions (sheet names must map)
    out_pdf: output path for the figure PDF (or stem when split=True)
    """
    if buffer_groups is None:
        buffer_groups = [
    {
        'group_key':'HEPES',
        'base':'HEP_base',
        'acid':'H_HEP_acid',         
        'display_name':'HEPES/H-HEPES',
        'factor':1e3,'unit':'mM',
        'hhpred_col_exact':'HEPES_A-',
        'hhpred_use_ph_align':True
    },
    {
        'group_key':'Acetate',
        'base':'Acetate_base',
        'acid':'Acetic_acid',         
        'display_name':'Acetate/Acetic acid',
        'factor':1e3,'unit':'mM',
        'hhpred_col_exact':'Mg_A-',
        'hhpred_use_ph_align':True
    },
    {
        'group_key':'Pi',
        'base':'Pi',
        'acid':'Pi',
        'display_name':r'$\mathrm{H-Pi}/\mathrm{H_{2}-Pi}$',
        'factor':1e3,'unit':'mM',
        'hhpred_col_exact':'Pi_A-',
        'hhpred_use_ph_align':True
    },
]

    matplotlib.rcParams['text.usetex'] = False

    xls = None
    hhpred_sheets = {}
    sheet_map = {}
    if hhpred_excel:
        xls = pd.ExcelFile(hhpred_excel)
        sheet_map = _sheet_map_from_excel(xls)

    def _sheet_for_label(label: str):
        return sheet_map.get(_normalize_name(label), None)

    def _plot_for_subset(subset_paths: dict[str, Path], title_suffix: str='ALL'):
        palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
        styles = {lbl: dict(color=palette[i%len(palette)], linestyle='-', lw=LINEWIDTH)
                  for i,lbl in enumerate(subset_paths)}

        n = len(buffer_groups)
        special = (n==3 and cols==2)

        if special:
            fig = plt.figure(figsize=figsize)
            gs = fig.add_gridspec(2, 4, height_ratios=[1, 1])
            axes = [fig.add_subplot(gs[0, 0:2]),
                    fig.add_subplot(gs[0, 2:4]),
                    fig.add_subplot(gs[1, 1:3])]
        else:
            rows = int(np.ceil(n/cols))
            fig, axes_arr = plt.subplots(rows, cols, figsize=figsize)
            axes = axes_arr.flatten()

        for idx,(ax,grp) in enumerate(zip(axes,buffer_groups)):
            fk = float(grp.get('factor',1.0))
            base_name = grp['base']; acid_name = grp['acid']
            display = grp.get('display_name', base_name)
            hh_col = grp.get('hhpred_col_exact')
            use_ph = bool(grp.get('hhpred_use_ph_align', False))
            key = grp.get('group_key','').lower()

            base_ls = '-'; acid_ls = ':'

            # panel label + title
            letter = chr(65 + idx)
            ax.text(0.0,  PANEL_TITLE_Y, f'({letter})', transform=ax.transAxes,
                    ha='left', va='bottom', fontsize=PANEL_LABEL_SIZE, fontweight='bold', clip_on=False)
            ax.text(TITLE_X_OFF, PANEL_TITLE_Y, f'{display}', transform=ax.transAxes,
                    ha='left', va='bottom', fontsize=PANEL_LABEL_SIZE, clip_on=False)

            for label, path in subset_paths.items():
                df = pd.read_csv(path)
                # time column
                tcol = pick_column(df, include_keys=('time','min')) or 'time_min'
                t = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
                # optional predicted pH column for alignment
                ph_col = pick_column(df, include_keys=('ph','pred')) or pick_column(df, include_keys=('pH_pred',)) or None
                ph_series = pd.to_numeric(df[ph_col], errors="coerce") if ph_col in df.columns else None

                # kinetic base
                if base_name:
                    bcol = pick_column(df, include_keys=(base_name.lower(),), exclude_keys=('minus','plus')) or base_name
                    yb = pd.to_numeric(df[bcol], errors="coerce").to_numpy() * fk
                    if smooth:
                        yb = pd.Series(yb).rolling(window=window, center=True, min_periods=1).mean().to_numpy()
                    yb = np.clip(yb, 0, None)
                    bs = styles[label].copy(); bs['linestyle'] = base_ls
                    m = np.isfinite(t) & np.isfinite(yb)
                    if m.any():
                        ax.plot(t[m], yb[m], label=f"{_short_label(label)} – Base", **bs)

                # kinetic acid
                ya = None
                if key=='acetate' and base_name:
                    bcol = pick_column(df, include_keys=(base_name.lower(),), exclude_keys=('minus','plus')) or base_name
                    base_vals = pd.to_numeric(df[bcol], errors="coerce").to_numpy()
                    ya = (0.042 - base_vals) * fk  # uses your acetate total assumption
                elif acid_name:
                    acol = pick_column(df, include_keys=(acid_name.lower(),), exclude_keys=('minus','plus')) or acid_name
                    ya = pd.to_numeric(df[acol], errors="coerce").to_numpy() * fk
                if ya is not None:
                    if smooth:
                        ya = pd.Series(ya).rolling(window=window, center=True, min_periods=1).mean().to_numpy()
                    ya = np.clip(ya, 0, None)
                    as_ = styles[label].copy(); as_['linestyle'] = acid_ls
                    m = np.isfinite(t) & np.isfinite(ya)
                    if m.any():
                        ax.plot(t[m], ya[m], label=f"{_short_label(label)} – Acid", **as_)

                if xls is not None and hh_col and key!='pi':
                    sheet = _sheet_for_label(label)
                    if sheet:
                        if sheet not in hhpred_sheets:
                            df_hh = pd.read_excel(xls, sheet_name=sheet)
                            df_hh.columns = [c.strip() for c in df_hh.columns]
                            hhpred_sheets[sheet] = df_hh
                        else:
                            df_hh = hhpred_sheets[sheet]

                        # pH alignment (optional)
                        if use_ph and ph_series is not None:
                            hh_ph = pick_column(df_hh, include_keys=('pH',)) or pick_column(df_hh, include_keys=('ph',))
                            if hh_ph:
                                hh_vals = pd.to_numeric(df_hh[hh_ph], errors="coerce").to_numpy()
                                tt = []
                                phv = ph_series.to_numpy(dtype=float)
                                for v in hh_vals:
                                    if not np.isfinite(v): 
                                        tt.append(np.nan); 
                                        continue
                                    diffs = np.abs(phv - v)
                                    if np.all(~np.isfinite(diffs)):
                                        tt.append(np.nan)
                                    else:
                                        tt.append(t[np.nanargmin(diffs)])
                                t_h = np.array(tt, dtype=float)
                            else:
                                t_h = None
                        else:
                            t_h = None

                        tgt = hh_col.strip().lower()
                        matched_col = None
                        for c in df_hh.columns:
                            if c.strip().lower() == tgt:
                                matched_col = c; break
                        if not matched_col:
                            norm_t = re.sub(r'[^a-z0-9]+','', tgt)
                            for c in df_hh.columns:
                                if norm_t in re.sub(r'[^a-z0-9]+','', c.lower()):
                                    matched_col = c; break

                        if matched_col:
                            hh_series = pd.to_numeric(df_hh[matched_col], errors="coerce").to_numpy()
                            if t_h is None or np.all(~np.isfinite(t_h)):
                                # fall back to uniform placement across [0, t_max]
                                t_max = np.nanmax(t) if np.any(np.isfinite(t)) else 120.0
                                t_h = np.linspace(0, t_max, len(hh_series))
                            L = min(len(t_h), len(hh_series))
                            tp = np.asarray(t_h[:L], dtype=float)
                            y_p = np.clip(hh_series[:L] * fk, 0, None)
                            color = styles[label]['color']
                            m = np.isfinite(tp) & np.isfinite(y_p)
                            if m.any():
                                ax.plot(tp[m], y_p[m], linestyle='none', marker='o', markersize=7,
                                        markeredgecolor=color, markerfacecolor='none')

                            if key in ('hepes','acetate'):
                                tot = 0.04 if key=='hepes' else 0.042
                                yc = np.clip((tot - hh_series[:L]) * fk, 0, None)
                                m2 = np.isfinite(tp) & np.isfinite(yc)
                                if m2.any():
                                    ax.plot(tp[m2], yc[m2], linestyle='none', marker='s', markersize=7,
                                            markeredgecolor=color, markerfacecolor='none')

            ax.set_ylim(bottom=0)
            _apply_axis_style(ax, grp.get("unit","mM"))

        if not special:
            for ax in axes[len(buffer_groups):]:
                ax.remove()

        handles, labels = [], []
        def _add_kinetic(lbl):
            if not lbl: return
            short = _short_label(lbl)
            color = styles[lbl]['color']
            handles.append(Line2D([0],[0], color=color, lw=LINEWIDTH, linestyle='-'))
            labels.append(f"Kinetic (UKF) Base – {short} mRNA")
            handles.append(Line2D([0],[0], color=color, lw=LINEWIDTH, linestyle=':'))
            labels.append(f"Kinetic (UKF) Acid – {short} mRNA")

        egfp_lbl = next((l for l in subset_paths if 'egfp' in l.lower()), None)
        csp_lbl  = next((l for l in subset_paths if 'csp'  in l.lower()), None)
        _add_kinetic(egfp_lbl); _add_kinetic(csp_lbl)

        if xls is not None:
            for lbl in (egfp_lbl, csp_lbl):
                if not lbl: continue
                color = styles[lbl]['color']; short = _short_label(lbl)
                handles += [
                    Line2D([0],[0], color=color, marker='o', linestyle='none', markersize=8,
                           markerfacecolor='none', markeredgecolor=color),
                    Line2D([0],[0], color=color, marker='s', linestyle='none', markersize=8,
                           markerfacecolor='none', markeredgecolor=color),
                ]
                labels += [f"H–H Prediction Base – {short} mRNA",
                           f"H–H Prediction Acid – {short} mRNA"]

        fig.legend(handles=handles, labels=labels,
                   loc='upper center', ncol=2, frameon=False, fontsize=LEGEND_FS,
                   markerscale=1.2, bbox_to_anchor=(0.5, LEGEND_Y),
                   borderaxespad=0.0, handletextpad=0.5,
                   columnspacing=1.2, labelspacing=0.3)

        fig.subplots_adjust(left=LEFT_FIG, right=RIGHT_FIG, bottom=BOT_FIG, top=TOP_FIG,
                            wspace=WSPACE, hspace=HSPACE)

        out_dir = Path(out_pdf).parent
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = Path(out_pdf)
        fig.savefig(out_path, dpi=600, bbox_inches='tight')
        plt.close(fig)

    if split_by_tris_hepes:
        hep = {k:v for k,v in dataset_paths.items() if 'HEPES' in k.upper()}
        tris= {k:v for k,v in dataset_paths.items() if 'TRIS'  in k.upper()}

        if hep:
            p = Path(out_pdf); out_hep = p.with_name(p.stem.replace("ALL","HEPES") + p.suffix)
            _plot_for_subset(hep,  title_suffix='HEPES')
            Path(out_pdf).rename(out_hep)
        if tris:
            p = Path(out_pdf); out_tris = p.with_name(p.stem.replace("ALL","TRIS") + p.suffix)
            _plot_for_subset(tris, title_suffix='TRIS')
            Path(out_pdf).rename(out_tris)
    else:
        _plot_for_subset(dataset_paths, title_suffix='ALL')


# ------------------- CLI -------------------
def main():
    ap = argparse.ArgumentParser(description="Buffer base/acid plots from UKF CSVs with optional HH overlays.")
    ap.add_argument("--in-dir", default="reports/softsensor_run", help="Directory with UKF CSVs.")
    ap.add_argument("--hh-excel", default=None, help="Path to HH workbook (e.g., data/H_H_NTP.xlsx).")
    ap.add_argument("--out", default="figures/buffers/buffers_ALL.pdf", help="Output PDF path.")
    ap.add_argument("--split", action="store_true", help="Split into HEPES/TRIS figures.")
    ap.add_argument("--smooth", action="store_true")
    ap.add_argument("--window", type=int, default=5)
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    dataset_paths = {
        "eGFP + HEPES": in_dir / CSV_DEFAULTS["eGFP + HEPES"],
        "eGFP + TRIS":  in_dir / CSV_DEFAULTS["eGFP + TRIS"],
        "CSP + HEPES":  in_dir / CSV_DEFAULTS["CSP + HEPES"],
        "CSP + TRIS":   in_dir / CSV_DEFAULTS["CSP + TRIS"],
    }

    plot_buffer_groups(
        dataset_paths=dataset_paths,
        smooth=args.smooth,
        window=args.window,
        cols=2,
        figsize=FIGSIZE,
        split_by_tris_hepes=args.split,
        hhpred_excel=args.hh_excel,
        out_pdf=args.out,
    )
    print(f"[OK] wrote {args.out}")

if __name__ == "__main__":
    main()
