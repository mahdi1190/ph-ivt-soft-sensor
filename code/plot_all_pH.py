#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_all_pH.py
==============

Four-panel figure (2×2) showing:
- pH traces (ph1/2/3 scatter + per-time mean line)
- RNA concentration from UV spectrophotometry (rna1/2/3) with error bars
- RNA concentration from Fluorometry (qubit1/2) with error bars
- Temperature (°C)

Source data: data/all_data_processed.xlsx with sheets:
    egfphepes, egfptris, csphepes, csptris
Required columns per sheet (any can have blanks/NaNs):
    time, ph1, ph2, ph3, temp, rna1, rna2, rna3, qubit1, qubit2

Outputs:
    figures/pH_quads/pH_rna_temp_quad.png
    figures/pH_quads/pH_rna_temp_quad.pdf

Author: Mahdi Ahmed (2025)
License: MIT
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Tuple, List, Dict
import argparse
import numpy as np
import pandas as pd
import matplotlib
if not matplotlib.get_backend().lower().startswith("qt"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams['axes.titlelocation'] = 'left'
COL = {
    'ph'  : '#003366',  # dark blue
    'rna' : '#006633',  # dark green (used for both UV & Fluorometry)
    'temp': '#CC6600',  # dark orange
}

TITLE_SIZE      = 12
AX_LABEL_SIZE   = 11
TICK_LABEL_SIZE = 10

FIG_W, FIG_H = 7.3, 6.5
TOP, BOT, LEFT, RIGHT = 0.86, 0.10, 0.10, 0.98
WSPACE, HSPACE = 0.30, 0.35
LEG_FONTSIZE = 10

NA_STRINGS = ["", " ", "NA", "N/A", "na", "NaN", "-", "--"]

RUNS: List[Tuple[str, str]] = [
    ("egfphepes", "eGFP + HEPES"),
    ("egfptris",  "eGFP + TRIS"),
    ("csphepes",  "CSP + HEPES"),
    ("csptris",   "CSP + TRIS"),
]

def _clean(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [c.strip().replace(" ", "").lower() for c in out.columns]
    # coerce known numeric columns if they exist
    for c in ["time","ph1","ph2","ph3","temp","rna1","rna2","rna3","qubit1","qubit2"]:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out

def _finite_mask(*arrays) -> np.ndarray:
    mask = np.ones_like(np.asarray(arrays[0], dtype=float), dtype=bool)
    for a in arrays:
        mask &= np.isfinite(np.asarray(a, dtype=float))
    return mask

def _pick_cols(df: pd.DataFrame, names: Iterable[str]) -> List[str]:
    return [c for c in names if c in df.columns]

# ------------------- panel plot -------------------
def plot_sheet(ax: plt.Axes, df: pd.DataFrame, title: str, unc_scale: float = 1.0) -> Tuple[List, List]:
    """
    Draw one panel:
      - pH scatter clouds (ph1/2/3) + mean line
      - UV RNA (rna1/2/3) with error bars
      - Fluorometry RNA (qubit1/2) with error bars
      - Temperature on tertiary axis
    Returns handles, labels for legend aggregation.
    """
    # time axis (drop rows without valid time)
    if "time" not in df.columns:
        raise ValueError("Sheet missing 'time' column.")
    m_t = np.isfinite(df["time"].to_numpy())
    df = df.loc[m_t].sort_values("time")
    t = df["time"].to_numpy()

    handles, labels = [], []

    # ---- pH ----
    ph_cols = _pick_cols(df, ["ph1","ph2","ph3"])
    if ph_cols:
        # light scatter for each available pH replicate
        for col in ph_cols:
            m = _finite_mask(t, df[col].to_numpy())
            if m.any():
                h = ax.scatter(t[m], df[col].to_numpy()[m], s=12, alpha=0.25, color=COL['ph'])
        ph_mean = df[ph_cols].mean(axis=1, skipna=True).to_numpy()
        m = _finite_mask(t, ph_mean)
        if m.any():
            h_line, = ax.plot(t[m], ph_mean[m], '-', lw=2, color=COL['ph'], label="pH")
            handles.append(h_line); labels.append("pH")
        ax.set_ylabel("pH", fontsize=AX_LABEL_SIZE, color=COL['ph'])
        ax.tick_params(axis='y', labelcolor=COL['ph'], labelsize=TICK_LABEL_SIZE)

    # ---- RNA (UV + Fluorometry) on twin axis ----
    ax2 = ax.twinx()
    # UV: rna1/2/3
    uv_cols = _pick_cols(df, ["rna1","rna2","rna3"])
    if uv_cols:
        uv_vals = df[uv_cols]
        m_uv = uv_vals.notna().any(axis=1) & np.isfinite(t)
        t_uv = t[m_uv.to_numpy()]
        if t_uv.size:
            uv_mean = uv_vals[m_uv].mean(axis=1, skipna=True).to_numpy()
            uv_min  = uv_vals[m_uv].min(axis=1, skipna=True).to_numpy()
            uv_max  = uv_vals[m_uv].max(axis=1, skipna=True).to_numpy()
            uv_err  = ((uv_max - uv_min) / 2.0) * float(unc_scale)
            h_uv = ax2.errorbar(t_uv, uv_mean, yerr=uv_err, fmt='-o', capsize=3, ms=5,
                                color=COL['rna'], label="UV spectrophotometry", alpha=0.95)
            handles.append(h_uv); labels.append("UV spectrophotometry")

    # Fluorometry: qubit1/2
    if "qubit1" in df.columns and "qubit2" in df.columns:
        qv = df[["qubit1","qubit2"]]
        m_q = qv.notna().any(axis=1) & np.isfinite(t)
        t_q = t[m_q.to_numpy()]
        if t_q.size:
            q_mean = qv[m_q].mean(axis=1, skipna=True).to_numpy()
            q_min  = qv[m_q].min(axis=1, skipna=True).to_numpy()
            q_max  = qv[m_q].max(axis=1, skipna=True).to_numpy()
            q_err  = ((q_max - q_min) / 2.0) * float(unc_scale)
            h_q = ax2.errorbar(t_q, q_mean, yerr=q_err, fmt='--d', capsize=3, ms=5,
                               color=COL['rna'], label="Fluorometry", alpha=0.85)
            handles.append(h_q); labels.append("Fluorometry")

    ax2.set_ylabel("RNA conc. (g/L)\n(UV & Fluorometry)", fontsize=AX_LABEL_SIZE, color=COL['rna'])
    ax2.tick_params(axis='y', labelcolor=COL['rna'], labelsize=TICK_LABEL_SIZE)

    # ---- Temperature on a third axis ----
    ax3 = ax.twinx()
    ax3.spines['right'].set_position(('outward', 48))
    if "temp" in df.columns:
        m_tmp = _finite_mask(t, df["temp"].to_numpy())
        if m_tmp.any():
            h_tmp, = ax3.plot(t[m_tmp], df["temp"].to_numpy()[m_tmp], '--', lw=2,
                              color=COL['temp'], label="Temperature", alpha=0.9)
            handles.append(h_tmp); labels.append("Temperature")
    ax3.set_ylabel("Temperature (°C)", fontsize=AX_LABEL_SIZE, color=COL['temp'])
    ax3.set_ylim(30, 40)
    ax3.tick_params(axis='y', labelcolor=COL['temp'], labelsize=TICK_LABEL_SIZE)

    ax.set_xlabel("Time (minutes)", fontsize=AX_LABEL_SIZE)
    ax.set_xlim(left=0)
    ax.tick_params(labelsize=TICK_LABEL_SIZE)
    ax2.tick_params(labelsize=TICK_LABEL_SIZE)
    ax3.tick_params(labelsize=TICK_LABEL_SIZE)

    return handles, labels

# ------------------- main quad -------------------
def quad_plot(excel_path: str | Path,
              sheets: List[str],
              titles: List[str],
              uncales: List[float] | None = None,
              output_dir: str | Path = "figures/pH_quads",
              output_name: str = "pH_rna_temp_quad") -> Dict[str, Path]:
    """Create the 2×2 figure for the given four sheets."""
    excel_path = Path(excel_path)
    out_dir = Path(output_dir); out_dir.mkdir(parents=True, exist_ok=True)
    if uncales is None:
        uncales = [1.0] * len(sheets)

    dfs = pd.read_excel(excel_path,
                        sheet_name=sheets,
                        na_values=NA_STRINGS,
                        keep_default_na=True)
    for k in dfs:
        dfs[k] = _clean(dfs[k])

    fig, axes = plt.subplots(2, 2, figsize=(FIG_W, FIG_H), sharex=True)
    fig.subplots_adjust(left=LEFT, right=RIGHT, top=TOP, bottom=BOT, wspace=WSPACE, hspace=HSPACE)

    all_h, all_l = [], []
    letter_offset = 0.012

    for idx, (ax, sheet, title, unc) in enumerate(zip(axes.flatten(), sheets, titles, uncales)):
        h, l = plot_sheet(ax, dfs[sheet], title, unc_scale=unc)
        all_h.extend(h); all_l.extend(l)

        letter = chr(65 + idx)  # 'A', 'B', ...
        pos = ax.get_position()
        fig.text(
            pos.x0, pos.y1 + letter_offset,
            rf'$(\mathbf{{{letter}}})$ {title}',
            transform=fig.transFigure,
            fontsize=TITLE_SIZE, fontweight='normal', va='bottom', ha='left'
        )

    order = ["pH", "UV spectrophotometry", "Fluorometry", "Temperature"]
    legend_handles, legend_labels = [], []
    for name in order:
        if name in all_l:
            i = all_l.index(name)
            legend_handles.append(all_h[i])
            legend_labels.append(name)

    fig.legend(legend_handles, legend_labels,
               loc="upper center", ncol=len(legend_labels), frameon=False,
               fontsize=LEG_FONTSIZE, bbox_to_anchor=(0.5, 0.96))

    png = out_dir / f"{output_name}.png"
    pdf = out_dir / f"{output_name}.pdf"
    fig.savefig(png, dpi=600, bbox_inches="tight")
    fig.savefig(pdf, dpi=600, bbox_inches="tight")
    plt.close(fig)
    return {"png": png, "pdf": pdf}

def main():
    ap = argparse.ArgumentParser(description="2×2 pH/RNA/Temp quad plot from all_data_processed.xlsx")
    ap.add_argument("--all-data", default="data/all_data_processed.xlsx",
                    help="Path to workbook with sheets egfphepes/egfptris/csphepes/csptris.")
    ap.add_argument("--out-dir", default="figures/pH_quads", help="Output directory.")
    ap.add_argument("--name", default="pH_rna_temp_quad", help="Output filename stem (no extension).")
    args = ap.parse_args()

    sheets = [s for s, _ in RUNS]
    titles = [t for _, t in RUNS]
    outs = quad_plot(args.all_data, sheets, titles, uncales=[1.0]*4,
                     output_dir=args.out_dir, output_name=args.name)
    print({k: str(v) for k, v in outs.items()})

if __name__ == "__main__":
    main()
