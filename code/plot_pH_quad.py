#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_ph_quad.py — 4× pH curves (experiment + UKF) from the Excel + reports CSVs.

Author: Mahdi Ahmed (2025)
License: MIT
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
if not matplotlib.get_backend().lower().startswith("qt"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# ---------- config ----------
ALL_DATA_XLSX = "data/all_data_processed.xlsx"
UKF_DIR       = "reports/softsensor_run"
OUT_DIR       = "figures/ph"

PH_COL_MAP = {
    "egfphepes": "ph2",
    "egfptris":  "ph1",
    "csphepes":  "ph2",
    "csptris":   "ph3",
}

UKF_FILES = {
    "CSP + HEPES":  "ivt_ukf_results_csp_HEPES.csv",
    "CSP + TRIS":   "ivt_ukf_results_csp_TRIS.csv",
    "eGFP + HEPES": "ivt_ukf_results_egfp_HEPES.csv",
    "eGFP + TRIS":  "ivt_ukf_results_egfp_TRIS.csv",
}

# Sheet names for experiment
SHEET_FOR_LABEL = {
    "CSP + HEPES":  "csphepes",
    "CSP + TRIS":   "csptris",
    "eGFP + HEPES": "egfphepes",
    "eGFP + TRIS":  "egfptris",
}

# Plot style
CB_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]  # Okabe–Ito
MARKERS   = ["o", "s", "^", "d"]
SCALE_CANVAS = 1.0
BASE_W, BASE_H = 7.3, 7.3
LEGEND_FONTSIZE = 10
LEG1_Y, LEG2_Y, LEG3_Y = 0.77, 0.735, 0.70
TOP_RECT = 0.66  # headroom for three stacked legend rows
LINEWIDTH = 2.0
MARKERSIZE = 5

# Smoothing for UKF pH
SMOOTH_PREDICTIONS = True
SMOOTHING_WINDOW   = 3

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
})

# ---------- helpers ----------
def _norm(s: str) -> str:
    return "".join(ch for ch in str(s).lower() if ch.isalnum())

def _pick(df: pd.DataFrame, *candidates: str) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    cols_norm = [_norm(c) for c in df.columns]
    for cand in candidates:
        key = _norm(cand)
        for i, cn in enumerate(cols_norm):
            if key and key in cn:
                return df.columns[i]
    return None

def apply_plain_ticks(ax: plt.Axes):
    fmt = ScalarFormatter(useOffset=False, useMathText=False)
    fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)

def load_experimental_ph(all_data_xlsx: Path, sheet: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_excel(all_data_xlsx, sheet_name=sheet)
    tcol = _pick(df, "time")
    if tcol is None:
        raise KeyError(f"Sheet '{sheet}' needs a 'time' column.")
    # Choose mapped pH column; if missing, fall back to first available among ph1/ph2/ph3
    ph_pref = PH_COL_MAP.get(sheet, None)
    if ph_pref and ph_pref in df.columns:
        pcol = ph_pref
    else:
        pcol = _pick(df, "ph1", "ph2", "ph3")
    if pcol is None:
        raise KeyError(f"Sheet '{sheet}' needs a pH column (ph1/ph2/ph3).")

    t = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
    y = pd.to_numeric(df[pcol], errors="coerce").to_numpy()
    good = np.isfinite(t) & np.isfinite(y)
    return t[good], y[good]

def load_ukf_ph(csv_path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not csv_path.exists():
        return np.array([]), np.array([])
    df = pd.read_csv(csv_path)
    tcol = _pick(df, "time_min", "time", "minute")
    ycol = _pick(df, "pH_pred", "ph_pred", "pH", "ph")
    if tcol is None or ycol is None:
        return np.array([]), np.array([])
    t = pd.to_numeric(df[tcol], errors="coerce").to_numpy()
    y = pd.to_numeric(df[ycol], errors="coerce").to_numpy()
    good = np.isfinite(t) & np.isfinite(y)
    return t[good], y[good]

# ---------- main plotting ----------
def plot_ph_quad(
    all_data_path: str | Path = ALL_DATA_XLSX,
    ukf_dir: str | Path = UKF_DIR,
    out_dir: str | Path = OUT_DIR,
    smooth_predictions: bool = SMOOTH_PREDICTIONS,
    smoothing_window: int = SMOOTHING_WINDOW,
):
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)

    labels = ["CSP + HEPES", "CSP + TRIS", "eGFP + HEPES", "eGFP + TRIS"]
    # load all experimental first (order fixes legend grouping)
    exp_series = []
    for lbl in labels:
        sheet = SHEET_FOR_LABEL[lbl]
        t_exp, y_exp = load_experimental_ph(Path(all_data_path), sheet)
        exp_series.append((t_exp, y_exp, lbl))

    # load UKF preds
    pred_series = []
    for lbl in labels:
        csv_path = Path(ukf_dir) / UKF_FILES[lbl]
        t_pred, y_pred = load_ukf_ph(csv_path)
        if smooth_predictions and y_pred.size:
            y_pred = (pd.Series(y_pred)
                      .rolling(window=smoothing_window, min_periods=1, center=True)
                      .mean().to_numpy())
        pred_series.append((t_pred, y_pred, lbl))

    # figure
    fig, ax = plt.subplots(
        facecolor="white",
        figsize=(BASE_W * SCALE_CANVAS, BASE_H * SCALE_CANVAS)
    )

    # 1) Experimental curves (solid, coloured, markers)
    for (t_e, y_e, lbl), color, marker in zip(exp_series, CB_COLORS, MARKERS):
        ax.plot(
            t_e, y_e,
            linestyle="-", color=color, linewidth=1.8, alpha=0.9,
            marker=marker, markersize=MARKERSIZE, markevery=max(1, len(t_e)//24),
            label=lbl
        )

    # 2) Predicted curves (dashed, matching colour; single legend entry)
    for i, ((t_p, y_p, lbl), color) in enumerate(zip(pred_series, CB_COLORS)):
        if y_p.size:
            t_e, y_e, _ = exp_series[i]
            if y_e.size and y_p.size:
                y_p = y_p.copy()
                y_p[0] = y_e[np.argmin(np.abs(t_e - t_e.min()))]
            ax.plot(
                t_p if t_p.size else np.linspace(0, 120, len(y_p)),
                y_p,
                linestyle="--", color=color, linewidth=2.0, alpha=1.0,
                label="Kinetic Model (UKF)" if i == 0 else "_nolegend_"
            )

    ax.set_xlabel("Time (minutes)")
    ax.set_ylabel("pH")
    ax.set_xlim(0, 120)
    ymins = []; ymaxs = []
    for t_e, y_e, _ in exp_series:
        if y_e.size: ymins.append(np.nanmin(y_e)); ymaxs.append(np.nanmax(y_e))
    for t_p, y_p, _ in pred_series:
        if y_p.size: ymins.append(np.nanmin(y_p)); ymaxs.append(np.nanmax(y_p))
    if ymins and ymaxs:
        ylo = min(ymins) - 0.1
        yhi = max(ymaxs) + 0.1
        ax.set_ylim(ylo, yhi)
    apply_plain_ticks(ax)
    ax.grid(False)
    handles, labels_all = ax.get_legend_handles_labels()
    h_by_l = {l: h for h, l in zip(handles, labels_all)}

    row1 = ["eGFP + HEPES", "eGFP + TRIS"]
    row2 = ["CSP + HEPES", "CSP + TRIS"]
    row3 = ["Kinetic Model (UKF)"]

    fig.legend(
        [h_by_l[l] for l in row1 if l in h_by_l],
        [l for l in row1 if l in h_by_l],
        loc="upper center", bbox_to_anchor=(0.5, LEG1_Y),
        ncol=2, frameon=False, prop={'size': LEGEND_FONTSIZE},
        handlelength=1.8, handletextpad=0.6, columnspacing=1.6, labelspacing=0.6,
    )
    fig.legend(
        [h_by_l[l] for l in row2 if l in h_by_l],
        [l for l in row2 if l in h_by_l],
        loc="upper center", bbox_to_anchor=(0.5, LEG2_Y),
        ncol=2, frameon=False, prop={'size': LEGEND_FONTSIZE},
        handlelength=1.8, handletextpad=0.6, columnspacing=1.6, labelspacing=0.6,
    )
    if row3[0] in h_by_l:
        fig.legend(
            [h_by_l[row3[0]]], [row3[0]],
            loc="upper center", bbox_to_anchor=(0.5, LEG3_Y),
            ncol=1, frameon=False, prop={'size': LEGEND_FONTSIZE},
            handlelength=1.8, handletextpad=0.6, columnspacing=1.6, labelspacing=0.6,
        )

    plt.tight_layout(rect=(0, 0, 1, TOP_RECT))

    # Save
    out_pdf = Path(out_dir) / "ph_4curves_with_predictions_cb.pdf"
    out_png = Path(out_dir) / "ph_4curves_with_predictions_cb.png"
    fig.savefig(out_pdf, dpi=600, bbox_inches="tight")
    fig.savefig(out_png, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] wrote:\n  {out_pdf}\n  {out_png}")

if __name__ == "__main__":
    plot_ph_quad(
        all_data_path=ALL_DATA_XLSX,
        ukf_dir=UKF_DIR,
        out_dir=OUT_DIR,
        smooth_predictions=SMOOTH_PREDICTIONS,
        smoothing_window=SMOOTHING_WINDOW,
    )
