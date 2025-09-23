#!/usr/bin/env python3
"""
plot_NTP_quad.py
================

Eight-panel figure (4 runs × [Model vs HH]):
- Experimental data: data/all_data_processed.xlsx with sheets:
    egfphepes, egfptris, csphepes, csptris
- pH selection rules: egfphepes→ph2, csptris→ph3, csphepes→ph2, egfptris→ph1
- Model CSVs read from reports/softsensor_run/
- HH workbook (Excel) (H_H_NTP.xlsx)

Usage (CLI):
  python -m code.plot_NTP_quad \
    --all-data data/all_data_processed.xlsx \
    --hh-excel data/H_H_NTP.xlsx \
    --out-root figures/quad_ntp

Or import:
  from code.plot_NTP_quad import eight_panel
  eight_panel("data/all_data_processed.xlsx", "data/H_H_NTP.xlsx")

Deps:
  numpy pandas matplotlib scikit-learn openpyxl

Author: Mahdi Ahmed (2025)
License: MIT
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Tuple
import argparse
import warnings

import numpy as np
import pandas as pd
import matplotlib
if not matplotlib.get_backend().lower().startswith("qt"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.metrics import r2_score, mean_squared_error

# ------------------- styling & constants -------------------
NTPC = {'ATP':'#0072B2','UTP':'#D55E00','GTP':'#009E73','CTP':'#CC79A7'}
MRK  = {'ATP':'o','UTP':'s','GTP':'^','CTP':'d'}

TITLE_SIZE        = 10
AX_LABEL_SIZE     = 10
TICK_LABEL_SIZE   = 10

BASE_W, BASE_H = 7.3, 10.0
TOP_FIG        = 0.86
BOT_FIG        = 0.08
LEFT_FIG       = 0.11
RIGHT_FIG      = 0.995
WSPACE, HSPACE = 0.30, 0.43

LEGEND_FONTSIZE = 12
LEG1_Y, LEG2_Y  = 0.987, 0.955

PANEL_LABEL_Y   = 1.03
TITLE_X_OFFSET  = 0.10

TABLE_COL_WIDTHS = (0.10, 0.20)
TABLE_HEIGHT     = 0.48
TABLE_FONTSIZE   = 6
DX_TABLE         = -0.03
DY_TABLE         = -0.12

NA_STRINGS = ["", " ", "NA", "N/A", "na", "NaN", "-", "--"]

RUNS: List[Tuple[str, str]] = [
    ("egfphepes", "eGFP + HEPES"),
    ("egfptris",  "eGFP + TRIS"),
    ("csphepes",  "CSP + HEPES"),
    ("csptris",   "CSP + TRIS"),
]

PH_PICK: Dict[str, str] = {
    "egfphepes": "ph2",
    "csptris":   "ph3",
    "csphepes":  "ph2",
    "egfptris":  "ph1",
}

MODEL_PATHS: Dict[str, str] = {
    "egfphepes": "reports/softsensor_run/ivt_ukf_results_egfp_HEPES.csv",
    "egfptris":  "reports/softsensor_run/ivt_ukf_results_egfp_TRIS.csv",
    "csphepes":  "reports/softsensor_run/ivt_ukf_results_csp_HEPES.csv",
    "csptris":   "reports/softsensor_run/ivt_ukf_results_csp_TRIS.csv",
}

# ------------------- NA / numeric safety -------------------
def _standardize_na(df: pd.DataFrame) -> pd.DataFrame:
    return df.replace(NA_STRINGS, np.nan)

def _coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    for c in cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")
    return df

def _finite_mask(*arrays) -> np.ndarray:
    mask = np.ones_like(np.asarray(arrays[0], dtype=float), dtype=bool)
    for a in arrays:
        mask &= np.isfinite(np.asarray(a, dtype=float))
    return mask

def _interp_clean(x, y, x_new):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x_new = np.asarray(x_new, dtype=float)
    m = _finite_mask(x, y)
    if m.sum() < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    xw, yw = x[m], y[m]
    order = np.argsort(xw)
    xw, yw = xw[order], yw[order]
    xw_unique, idx = np.unique(xw, return_index=True)
    yw = yw[idx]
    return np.interp(x_new, xw_unique, yw, left=np.nan, right=np.nan)

# ------------------- experimental loader (one workbook) -------------------
def _clean_cols(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out.columns = [c.strip().replace(" ", "_") for c in out.columns]
    return out

def load_exp(all_data_path: str | Path, run_id: str) -> pd.DataFrame:
    """
    Standardise experimental data for one run/sheet.
    Returns columns: time_min, pH, RNA, ATP_tot, GTP_tot, CTP_tot, UTP_tot [, temp]
    """
    path = Path(all_data_path)
    run = run_id.lower().strip()
    if run not in PH_PICK:
        raise ValueError(f"Unknown run_id '{run_id}'. Expected one of {list(PH_PICK)}")
    ph_col = PH_PICK[run]

    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, sheet_name=run, na_values=NA_STRINGS, keep_default_na=True)
    else:
        # CSV fallback expects a run_id column
        df_all = pd.read_csv(path, na_values=NA_STRINGS, keep_default_na=True)
        if "run_id" not in df_all.columns:
            raise ValueError("CSV must include a 'run_id' column.")
        df = df_all[df_all["run_id"].astype(str).str.lower() == run].copy()

    df = _clean_cols(df)
    req = ["time","rna2","ATP_tot","GTP_tot","CTP_tot","UTP_tot", ph_col]
    missing = [c for c in req if c not in df.columns]
    if missing:
        raise ValueError(f"[{run}] missing columns: {missing}. Have: {list(df.columns)}")

    nums = ["time","rna2","ATP_tot","GTP_tot","CTP_tot","UTP_tot", ph_col]
    if "temp" in df.columns: nums.append("temp")
    df = _coerce_numeric(df, nums)

    before = len(df)
    df = df[np.isfinite(df["time"])]
    if before - len(df):
        warnings.warn(f"[{run}] dropped {before - len(df)} rows with invalid time")

    out = pd.DataFrame({
        "time_min": df["time"].astype(float),
        "pH":       df[ph_col].astype(float),
        "RNA":      df["rna2"].astype(float),
        "ATP_tot":  df["ATP_tot"].astype(float),
        "GTP_tot":  df["GTP_tot"].astype(float),
        "CTP_tot":  df["CTP_tot"].astype(float),
        "UTP_tot":  df["UTP_tot"].astype(float),
    })
    if "temp" in df.columns:
        out["temp"] = df["temp"].astype(float)

    # keep rows that have at least something measured (avoid pure-NaN lines)
    keep = _finite_mask(out["time_min"]) & (
        _finite_mask(out["ATP_tot"]) | _finite_mask(out["GTP_tot"]) |
        _finite_mask(out["CTP_tot"]) | _finite_mask(out["UTP_tot"]) |
        _finite_mask(out["pH"]) | _finite_mask(out["RNA"])
    )
    return out[keep].sort_values("time_min").reset_index(drop=True)

# ------------------- metrics & table -------------------
def _metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float,float]:
    return r2_score(y_true, y_pred), float(np.sqrt(mean_squared_error(y_true, y_pred)))

def _choose_table_bbox(ax, col_widths=TABLE_COL_WIDTHS, table_height=TABLE_HEIGHT):
    w = sum(col_widths); h = table_height
    fixed = {'tl': (0.06, 0.62), 'tr': (0.72, 0.62), 'bl': (0.06, 0.05), 'br': (0.72, 0.05)}
    candidates = [fixed['tl'], fixed['tr'], fixed['bl'], fixed['br']]
    xlim = ax.get_xlim(); ylim = ax.get_ylim()
    xr = (xlim[1] - xlim[0]) or 1.0; yr = (ylim[1] - ylim[0]) or 1.0

    xs, ys = [], []
    for ln in ax.get_lines():
        x = np.asarray(ln.get_xdata()); y = np.asarray(ln.get_ydata())
        if len(x) and len(y):
            xs.append((x - xlim[0]) / xr); ys.append((y - ylim[0]) / yr)
    if not xs:
        x0, y0 = candidates[0]; return (x0, y0, w, h)

    X = np.concatenate(xs); Y = np.concatenate(ys)
    def score(x0, y0):
        return np.count_nonzero((X >= x0) & (X <= x0 + w) & (Y >= y0) & (Y <= y0 + h))
    x0, y0 = min(candidates, key=lambda c: score(*c))
    return (x0, y0, w, h)

def _draw_table(ax, metrics):
    x0, y0, w, h = _choose_table_bbox(ax)
    x0 += DX_TABLE; y0 += DY_TABLE
    cell_text = [[f"{r2:.3f}", f"{rmse:.3f}"] for _, r2, rmse in metrics]
    avg_r2   = float(np.mean([r for _, r, _ in metrics])) if metrics else float("nan")
    avg_rmse = float(np.mean([e for _, _, e in metrics])) if metrics else float("nan")
    cell_text.append([f"{avg_r2:.3f}", f"{avg_rmse:.3f}"])
    row_labels = [ntp for ntp, _, _ in metrics] + ['Avg']
    col_labels = ['R²', 'RMSE (mM)']
    tbl = ax.table(cellText=cell_text, rowLabels=row_labels, colLabels=col_labels,
                   cellLoc='center', colWidths=list(TABLE_COL_WIDTHS),
                   bbox=(x0, y0, w, h))
    tbl.scale(1, 1.35)
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(TABLE_FONTSIZE)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor('black'); cell.set_linewidth(0.6)
        if 1 <= r <= len(metrics):
            cell.get_text().set_color(NTPC[row_labels[r-1]])
        if r == len(metrics) + 1:
            cell.get_text().set_fontweight('bold')
    return tbl

# ------------------- panels -------------------
def _time_col(df: pd.DataFrame) -> str:
    for k in df.columns:
        if "time" in k.lower() or "minute" in k.lower():
            return k
    raise ValueError("No time/minute column found.")

def plot_model_vs_exp(ax, model_csv: str | Path, exp_df: pd.DataFrame,
                      ntp_keys=("ATP_tot","UTP_tot","GTP_tot","CTP_tot"),
                      annotate_metrics=True, smooth_window=5):
    dfm = pd.read_csv(model_csv)
    t_mod = pd.to_numeric(dfm[_time_col(dfm)], errors="coerce").to_numpy()
    t_exp = exp_df["time_min"].to_numpy()

    lw, msize = 2, 6
    metrics = []
    for key in ntp_keys:
        ntp = key.split('_')[0].upper()
        if key not in dfm.columns or key not in exp_df.columns:
            continue
        y_mod = pd.to_numeric(dfm[key], errors="coerce").to_numpy()
        y_exp = pd.to_numeric(exp_df[key], errors="coerce").to_numpy()

        if smooth_window and smooth_window > 1:
            y_mod = pd.Series(y_mod).rolling(window=smooth_window, center=True, min_periods=1).mean().to_numpy()

        m_mod = _finite_mask(t_mod, y_mod)
        m_exp = _finite_mask(t_exp, y_exp)

        if m_mod.sum() >= 2:
            ax.plot(t_mod[m_mod], y_mod[m_mod], '-', lw=lw, color=NTPC[ntp])
        if m_exp.any():
            ax.plot(t_exp[m_exp], y_exp[m_exp], '--', lw=lw, marker=MRK[ntp], ms=msize, color=NTPC[ntp], alpha=0.75)

        if m_mod.sum() >= 2 and m_exp.sum() >= 2:
            y_mod_at_exp = _interp_clean(t_mod[m_mod], y_mod[m_mod], t_exp[m_exp])
            mpair = _finite_mask(y_mod_at_exp, y_exp[m_exp])
            if mpair.sum() >= 2:
                r2, rmse = _metrics(y_exp[m_exp][mpair], y_mod_at_exp[mpair])
                metrics.append((ntp, r2, rmse))

    if annotate_metrics and metrics:
        _draw_table(ax, metrics)

    ax.set_xlabel('Time (minutes)', fontsize=AX_LABEL_SIZE, labelpad=2)
    ax.set_ylabel('NTP Concentration (mM)', fontsize=AX_LABEL_SIZE, labelpad=2)
    ax.set_xlim(0, 120); ax.set_xticks(np.arange(0, 121, 20))
    ax.tick_params(labelsize=TICK_LABEL_SIZE, which='both', direction='out',
                   bottom=True, top=False, left=True, right=False,
                   labelbottom=True, labelleft=True)
    for s in ax.spines.values(): s.set_visible(True)

def plot_hh_vs_exp(ax, hh_excel_path: str | Path | None, hh_sheet: str,
                   exp_df: pd.DataFrame, ntp_keys=("ATP_tot","UTP_tot","GTP_tot","CTP_tot"),
                   annotate_metrics=True, smooth_window=5):
    if hh_excel_path is None:
        # echo experimental curves if HH not provided
        t_exp = exp_df["time_min"].to_numpy()
        for key in ntp_keys:
            ntp = key.split('_')[0].upper()
            y_exp = pd.to_numeric(exp_df[key], errors="coerce").to_numpy()
            m_exp = _finite_mask(t_exp, y_exp)
            if m_exp.any():
                ax.plot(t_exp[m_exp], y_exp[m_exp], '--', lw=2, marker=MRK[ntp], ms=6, color=NTPC[ntp], alpha=0.75)
        ax.set_xlabel('Time (minutes)'); ax.set_ylabel('NTP Concentration (mM)')
        ax.set_xlim(0, 120); ax.set_xticks(np.arange(0,121,20))
        return

    raw = pd.read_excel(hh_excel_path, sheet_name=hh_sheet, na_values=NA_STRINGS, keep_default_na=True)
    raw.columns = [c.strip() for c in raw.columns]
    t_col = next((c for c in raw.columns if "time" in c.lower()), None) or "Time"
    dfm = pd.DataFrame({"time": pd.to_numeric(raw[t_col], errors="coerce").to_numpy() if t_col in raw.columns else np.arange(len(raw))})

    remap = {"ATP_tot":"A_remaining_mM","UTP_tot":"U_remaining_mM","GTP_tot":"G_remaining_mM","CTP_tot":"C_remaining_mM"}
    for k, src in remap.items():
        if src not in raw.columns:
            raise ValueError(f"[HH] expected '{src}' in sheet '{hh_sheet}'. Have: {list(raw.columns)}")
        dfm[k] = pd.to_numeric(raw[src], errors="coerce").to_numpy()

    t_mod = dfm["time"].to_numpy()
    t_exp = exp_df["time_min"].to_numpy()

    lw, msize = 2, 6
    metrics = []
    for key in ntp_keys:
        ntp = key.split('_')[0].upper()
        y_mod = pd.to_numeric(dfm[key], errors="coerce").to_numpy()
        y_exp = pd.to_numeric(exp_df[key], errors="coerce").to_numpy()

        if smooth_window and smooth_window > 1:
            y_mod = pd.Series(y_mod).rolling(window=smooth_window, center=True, min_periods=1).mean().to_numpy()

        m_mod = _finite_mask(t_mod, y_mod)
        m_exp = _finite_mask(t_exp, y_exp)

        if m_mod.sum() >= 2:
            ax.plot(t_mod[m_mod], y_mod[m_mod], '-',  lw=lw, color=NTPC[ntp])
        if m_exp.any():
            ax.plot(t_exp[m_exp], y_exp[m_exp], '--', lw=lw, marker=MRK[ntp], ms=msize, color=NTPC[ntp], alpha=0.75)

        if m_mod.sum() >= 2 and m_exp.sum() >= 2:
            y_mod_at_exp = _interp_clean(t_mod[m_mod], y_mod[m_mod], t_exp[m_exp])
            mpair = _finite_mask(y_mod_at_exp, y_exp[m_exp])
            if mpair.sum() >= 2:
                r2, rmse = _metrics(y_exp[m_exp][mpair], y_mod_at_exp[mpair])
                metrics.append((ntp, r2, rmse))

    if annotate_metrics and metrics:
        _draw_table(ax, metrics)

    ax.set_xlabel('Time (minutes)', fontsize=AX_LABEL_SIZE, labelpad=2)
    ax.set_ylabel('NTP Concentration (mM)', fontsize=AX_LABEL_SIZE, labelpad=2)
    ax.set_xlim(0, 120); ax.set_xticks(np.arange(0, 121, 20))
    ax.tick_params(labelsize=TICK_LABEL_SIZE, which='both', direction='out',
                   bottom=True, top=False, left=True, right=False,
                   labelbottom=True, labelleft=True)
    for s in ax.spines.values(): s.set_visible(True)

# ------------------- main figure -------------------
def eight_panel(all_data: str | Path,
                hh_excel: str | Path | None = None,
                out_root: str | Path = "figures/quad_ntp",
                fig_width: float = BASE_W,
                fig_height: float = BASE_H,
                annotate_metrics: bool = True,
                smooth_window: int = 5,
                save_pdf: bool = True) -> Dict[str, Path]:
    all_data = Path(all_data)
    out_root = Path(out_root); out_root.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(4, 2, figsize=(fig_width, fig_height), sharex=False, sharey=True)
    axes = np.asarray(axes)
    left_letters  = ['A','B','C','D']
    right_letters = ['E','F','G','H']

    for r, (run_id, label) in enumerate(RUNS):
        exp_df = load_exp(all_data, run_id)
        model_csv = MODEL_PATHS[run_id]
        hh_sheet  = run_id  # sheet names match run ids

        axL = axes[r, 0]
        plot_model_vs_exp(axL, model_csv, exp_df,
                          annotate_metrics=annotate_metrics, smooth_window=smooth_window)
        axL.text(0.0, PANEL_LABEL_Y, f'({left_letters[r]})', transform=axL.transAxes,
                 fontsize=TITLE_SIZE, fontweight='bold', ha='left', va='bottom')
        axL.text(TITLE_X_OFFSET, PANEL_LABEL_Y, label, transform=axL.transAxes,
                 fontsize=TITLE_SIZE, ha='left', va='bottom')

        axR = axes[r, 1]
        plot_hh_vs_exp(axR, hh_excel, hh_sheet, exp_df,
                       annotate_metrics=annotate_metrics, smooth_window=smooth_window)
        axR.text(0.0, PANEL_LABEL_Y, f'({right_letters[r]})', transform=axR.transAxes,
                 fontsize=TITLE_SIZE, fontweight='bold', ha='left', va='bottom')
        axR.text(TITLE_X_OFFSET, PANEL_LABEL_Y, label, transform=axR.transAxes,
                 fontsize=TITLE_SIZE, ha='left', va='bottom')

    # figure legends
    style_handles = [
        Line2D([], [], linestyle='-',  linewidth=2, color='0.3', label='Model (solid)'),
        Line2D([], [], linestyle='--', linewidth=2, color='0.3', label='Experimental (dashed)'),
    ]
    ntp_handles = [Line2D([], [], marker=MRK[k], linestyle='None', markersize=7,
                          label=k, color=NTPC[k]) for k in ['ATP','UTP','GTP','CTP']]

    fig.legend(handles=style_handles,
               bbox_to_anchor=(0.5, LEG1_Y), loc='upper center',
               ncol=2, frameon=False, fontsize=LEGEND_FONTSIZE,
               handlelength=1.6, handletextpad=0.6, columnspacing=1.2, labelspacing=0.4)
    fig.legend(handles=ntp_handles,
               bbox_to_anchor=(0.5, LEG2_Y), loc='upper center',
               ncol=4, frameon=False, fontsize=LEGEND_FONTSIZE,
               handlelength=0.0, handletextpad=0.6, columnspacing=1.2, labelspacing=0.4)

    fig.subplots_adjust(left=LEFT_FIG, right=RIGHT_FIG, bottom=BOT_FIG, top=TOP_FIG,
                        wspace=WSPACE, hspace=HSPACE)

    # column headers
    left_col  = [axes[i, 0] for i in range(axes.shape[0])]
    right_col = [axes[i, 1] for i in range(axes.shape[0])]
    x_left   = min(ax.get_position().x0 for ax in left_col)
    x_right  = min(ax.get_position().x0 for ax in right_col)
    grid_top = max(ax.get_position().y1 for ax in axes.ravel())
    y_header = min(0.995, grid_top + 0.023)
    fig.text(x_left,  y_header, "Kinetic Model (UKF)", ha='left', va='bottom',
             fontsize=TITLE_SIZE, fontweight='bold')
    fig.text(x_right, y_header, "Henderson–Hasselbalch", ha='left', va='bottom',
             fontsize=TITLE_SIZE, fontweight='bold')

    png = out_root / 'quad_ntp_profiles_8panel.png'
    pdf = out_root / 'quad_ntp_profiles_8panel.pdf'
    fig.savefig(png, dpi=600, bbox_inches='tight')
    if save_pdf:
        fig.savefig(pdf, dpi=600, bbox_inches='tight')
    plt.close(fig)
    return {"png": png, "pdf": pdf if save_pdf else None}

# ------------------- CLI -------------------
def main():
    p = argparse.ArgumentParser(description="8-panel NTP plots (Model vs HH) from single experimental workbook")
    p.add_argument("--all-data", required=True, help="Path to data/all_data_processed.xlsx (or CSV with run_id).")
    p.add_argument("--hh-excel", default="data/H_H_NTP.xlsx", help="Path to HH workbook (set 'None' to skip HH).")
    p.add_argument("--out-root", default="figures/quad_ntp", help="Output folder.")
    p.add_argument("--height", type=float, default=BASE_H, help="Figure height; width fixed at 7.3 in.")
    p.add_argument("--no-metrics", action="store_true")
    p.add_argument("--smooth-window", type=int, default=5)
    p.add_argument("--no-pdf", action="store_true")
    args = p.parse_args()

    hh = None if (args.hh_excel.lower() == "none") else args.hh_excel
    outs = eight_panel(
        all_data=args.all_data,
        hh_excel=hh,
        out_root=args.out_root,
        fig_height=float(args.height),
        annotate_metrics=not args.no_metrics,
        smooth_window=int(args.smooth_window),
        save_pdf=not args.no_pdf,
    )
    print({k: str(v) if v is not None else None for k, v in outs.items()})

if __name__ == "__main__":
    eight_panel("data/all_data_processed.xlsx", "data/H_H_NTP.xlsx")
