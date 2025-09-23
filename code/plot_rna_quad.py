#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_rna_quad.py 
Reads:
  • Experimental RNA from  data/all_data_processed.xlsx (sheets: egfphepes, egfptris, csphepes, csptris), column 'rna2'
  • UKF model CSVs from     reports/softsensor_run/
  • HH overlay from         data/HH.xlsx  (legacy: single sheet with columns egfp/egfp2/csp/csp2)
                             or data/H_H_NTP.xlsx (fallback: multi-sheet; any sheet with RNA-like column)

Outputs: figures/rna/quad_rna_profiles.(pdf|png)
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

# ---------- styling ----------
COL = {'exp':'#D55E00','ukf':'#0072B2','hh':'#009E73','unc':'#999999'}
AXFS=10; TICKFS=10; LABFS=10
LWS=2.2; MS=7
PANEL_Y=1.025

# canvas / layout
BASE_W, BASE_H = 7.3, 5.0
TOP=0.82; WSPACE=0.36; HSPACE=0.46
LEG_Y=0.985

# file structure defaults
ALL_DATA_XLSX = "data/all_data_processed.xlsx"
UKF_DIR       = "reports/softsensor_run"
OUT_DIR       = "figures/rna"

# mapping (pretty label -> sheet name / csv name)
SHEET_MAP = {
    "eGFP + HEPES": "egfphepes",
    "eGFP + TRIS":  "egfptris",
    "CSP + HEPES":  "csphepes",
    "CSP + TRIS":   "csptris",
}
CSV_DEFAULTS = {
    "eGFP + HEPES": "ivt_ukf_results_egfp_HEPES.csv",
    "eGFP + TRIS":  "ivt_ukf_results_egfp_TRIS.csv",
    "CSP + HEPES":  "ivt_ukf_results_csp_HEPES.csv",
    "CSP + TRIS":   "ivt_ukf_results_csp_TRIS.csv",
}
RNA_COL = "rna2"  # use rna2 for all sheets

# accept both legacy and new HH sheet namings
_HH_SHEET_ALIASES = {
    "egfphepes": "egfp",
    "egfptris" : "egfp2",
    "csphepes" : "csp",
    "csptris"  : "csp2",
}

def _norm(s: str) -> str:
    """normalize for matching: lowercase, alnum only."""
    return "".join(ch for ch in str(s).lower() if ch.isalnum())

def _pick(df: pd.DataFrame, *candidates: str) -> str | None:
    """Return the first exact (or loose) match from candidates present in df.columns."""
    for c in candidates:
        if c in df.columns: return c
    # loose (normalized contains)
    cols_norm = [_norm(c) for c in df.columns]
    for cand in candidates:
        k = _norm(cand)
        for i, cn in enumerate(cols_norm):
            if k and k in cn:
                return df.columns[i]
    return None

def apply_plain_ticks(ax: plt.Axes):
    """Force plain numeric y ticks (no scientific 1e±k or offsets)."""
    fmt = ScalarFormatter(useOffset=False, useMathText=False)
    fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)

def resolve_hh_path(candidates=("data/HH.xlsx", "data/H_H_NTP.xlsx")) -> Path | None:
    """Return the first HH workbook that exists, or None if neither is present."""
    for c in candidates:
        p = Path(c)
        if p.exists():
            return p
    return None

def load_exp_rna(all_data_xlsx: Path, sheet: str) -> tuple[np.ndarray, np.ndarray]:
    """Load experimental RNA (g/L) from the unified Excel. Uses 'rna2' and 'time'."""
    df = pd.read_excel(all_data_xlsx, sheet_name=sheet)
    tcol = _pick(df, "time")
    rcol = _pick(df, RNA_COL)
    if tcol is None or rcol is None:
        raise KeyError(f"Sheet '{sheet}': need 'time' and '{RNA_COL}' columns.")
    m = df[rcol].notna()
    t = pd.to_numeric(df.loc[m, tcol], errors="coerce").to_numpy()
    y = pd.to_numeric(df.loc[m, rcol], errors="coerce").to_numpy()
    good = np.isfinite(t) & np.isfinite(y)
    return t[good], y[good]

def load_hh_series(hh_path: Path | None, sheet_key: str) -> np.ndarray | None:
    """
    Robust HH loader:
      1) LEGACY (original): first sheet with columns 'egfp','egfp2','csp','csp2'
         -> pick the column matching alias for 'sheet_key'.
      2) FALLBACK: multi-sheet workbook -> find a sheet matching alias or sheet_key,
         then pick RNA-like column: ('RNA','rna','rna_hh','RNA_pred','pred') or first numeric column.
    Returns 1D array or None.
    """
    if not hh_path or not Path(hh_path).exists():
        return None

    alias = _HH_SHEET_ALIASES.get(sheet_key.lower(), sheet_key)
    alias_norm = _norm(alias)

    try:
        # (1) Try legacy: read FIRST sheet only and look for columns
        df0 = pd.read_excel(hh_path, sheet_name=0)
        # direct column match
        if alias in df0.columns:
            s = pd.to_numeric(df0[alias], errors="coerce").to_numpy()
            s = s[np.isfinite(s)]
            if s.size: return s
        # normalized fuzzy column match
        for col in df0.columns:
            if _norm(col) == alias_norm:
                s = pd.to_numeric(df0[col], errors="coerce").to_numpy()
                s = s[np.isfinite(s)]
                if s.size: return s
    except Exception:
        pass  # if any issue, fall through to multi-sheet scan

    # (2) Fallback: scan sheets
    xls = pd.ExcelFile(hh_path)
    # choose sheet by alias or original key (exact or fuzzy)
    chosen = None
    for name in xls.sheet_names:
        if _norm(name) == alias_norm or _norm(name) == _norm(sheet_key):
            chosen = name; break
    if chosen is None:
        # fuzzy contains
        for name in xls.sheet_names:
            if alias_norm in _norm(name) or _norm(sheet_key) in _norm(name):
                chosen = name; break
    if chosen is None:
        return None

    df = pd.read_excel(hh_path, sheet_name=chosen)
    # preferred RNA-like columns
    for cand in ("RNA","rna","rna_hh","RNA_pred","pred", alias):
        c = _pick(df, cand)
        if c:
            s = pd.to_numeric(df[c], errors="coerce").to_numpy()
            s = s[np.isfinite(s)]
            if s.size: return s
    # fallback: first numeric column
    for col in df.columns:
        s = pd.to_numeric(df[col], errors="coerce").to_numpy()
        s = s[np.isfinite(s)]
        if s.size: return s
    return None

def plot_panel(ax: plt.Axes,
               ukf_csv: Path,
               t_exp: np.ndarray, y_exp: np.ndarray,
               hh_vals: np.ndarray | None,
               letter: str,
               show_uncertainty: bool,
               unc_scale: float | None,
               smooth: bool,
               window: int,
               plot_model_at_exp: bool):
    """Draw one panel with experiment, UKF, and (optional) HH + CI."""
    dfm = pd.read_csv(ukf_csv)
    tcol = _pick(dfm, "time_min", "time", "minute")
    ycol = _pick(dfm, "RNA")
    if tcol is None or ycol is None:
        raise KeyError(f"{ukf_csv}: need time_min/time and RNA columns.")

    t_mod = pd.to_numeric(dfm[tcol], errors="coerce").to_numpy()
    y_mod = pd.to_numeric(dfm[ycol], errors="coerce").to_numpy()
    good = np.isfinite(t_mod) & np.isfinite(y_mod)
    t_mod = t_mod[good]; y_mod = y_mod[good]

    if smooth:
        y_mod = pd.Series(y_mod).rolling(window=window, center=True, min_periods=1).mean().to_numpy()

    # Uncertainty band (off by default)
    if show_uncertainty:
        lo = _pick(dfm, "RNA_minus_2σ","RNA_minus_2Ïσ","RNA_minus_2s","RNA_minus")
        hi = _pick(dfm, "RNA_plus_2σ","RNA_plus_2Ïσ","RNA_plus_2s","RNA_plus")
        if lo and hi:
            ylo = pd.to_numeric(dfm[lo], errors="coerce").to_numpy()[good]
            yhi = pd.to_numeric(dfm[hi], errors="coerce").to_numpy()[good]
            if smooth:
                ylo = pd.Series(ylo).rolling(window=window, center=True, min_periods=1).mean().to_numpy()
                yhi = pd.Series(yhi).rolling(window=window, center=True, min_periods=1).mean().to_numpy()
            if unc_scale is None: unc_scale = 1.0
            band_lo = np.clip(y_mod - unc_scale*(y_mod - ylo), 0, None)
            band_hi = np.clip(y_mod + unc_scale*(yhi - y_mod), 0, None)
            ax.fill_between(t_mod, band_lo, band_hi, color=COL['unc'], alpha=0.15, label="Uncertainty (±2σ)")

    # model curve
    y_mod_at_exp = np.interp(t_exp, t_mod, y_mod) if t_exp.size else np.array([])
    if plot_model_at_exp and t_exp.size:
        ax.plot(t_exp, y_mod_at_exp, color=COL['ukf'], lw=LWS, label="Kinetic Model (UKF)")
    else:
        ax.plot(t_mod, y_mod,     color=COL['ukf'], lw=LWS, label="Kinetic Model (UKF)")

    # HH overlay (optional)
    if hh_vals is not None and hh_vals.size:
        t_hh = np.linspace(t_mod.min() if t_mod.size else (t_exp.min() if t_exp.size else 0.0),
                           t_mod.max() if t_mod.size else (t_exp.max() if t_exp.size else 120.0),
                           len(hh_vals))
        if t_exp.size:
            y_hh_at_exp = np.interp(t_exp, t_hh, hh_vals)
            ax.plot(t_exp, y_hh_at_exp, color=COL['hh'], lw=LWS, linestyle="-.", label="Henderson–Hasselbalch")
        else:
            ax.plot(t_hh, hh_vals, color=COL['hh'], lw=LWS, linestyle="-.", label="Henderson–Hasselbalch")

    # experiment
    ax.plot(t_exp, y_exp, color=COL['exp'], lw=LWS, linestyle=":", marker="o", markersize=MS, label="Experiment")

    # metrics table (R², RMSE) comparing exp vs UKF & HH (at exp times)
    if t_exp.size:
        # UKF
        ss_res_mod = np.sum((y_exp - y_mod_at_exp)**2)
        ss_tot     = np.sum((y_exp - np.mean(y_exp))**2) if len(y_exp) > 1 else 1.0
        r2_mod     = 1 - ss_res_mod/ss_tot
        rmse_mod   = np.sqrt(ss_res_mod / max(len(y_exp), 1))
        # HH
        if hh_vals is not None and hh_vals.size:
            y_hh = np.interp(t_exp, t_hh, hh_vals)
            ss_res_hh = np.sum((y_exp - y_hh)**2)
            r2_hh     = 1 - ss_res_hh/ss_tot
            rmse_hh   = np.sqrt(ss_res_hh / max(len(y_exp), 1))
        else:
            r2_hh, rmse_hh = np.nan, np.nan

        cell_text = [[f'{r2_mod:.3f}', f'{rmse_mod:.3f}'],
                     [f'{(r2_hh if np.isfinite(r2_hh) else float("nan")):.3f}',
                      f'{(rmse_hh if np.isfinite(rmse_hh) else float("nan")):.3f}']]
        row_labels = ['UKF', 'H–H']
        col_labels = ['R²', 'RMSE (g/L)']
        tbl = ax.table(cellText=cell_text, rowLabels=row_labels, colLabels=col_labels,
                       cellLoc='center', bbox=(0.52, 0.12, 0.46, 0.28))
        tbl.auto_set_font_size(False); tbl.set_fontsize(8)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor('black')

    # axes + cosmetics
    ax.set_xlabel("Time (minutes)", fontsize=AXFS)
    ax.set_ylabel("RNA Yield (g/L)", fontsize=AXFS)
    ax.set_xlim(0, 120); ax.set_xticks(np.arange(0, 121, 20))
    ax.tick_params(labelsize=TICKFS)
    apply_plain_ticks(ax)  # plain numbers only
    ax.text(0.0, PANEL_Y, f"({letter})", transform=ax.transAxes,
            fontsize=LABFS, fontweight='bold', ha='left', va='bottom')

def quad_plot(
    in_dir: str | Path = UKF_DIR,
    all_data_path: str | Path = ALL_DATA_XLSX,
    hh_path: str | Path | None = None,
    out_dir: str | Path = OUT_DIR,
    show_uncertainty: bool = False,
    unc: float | None = 1.0,
    smooth: bool = True,
    window: int = 20,
    plot_model_at_exp: bool = False,
):
    """Build the 2×2 RNA figure and save to figures/rna/."""
    in_dir, out_dir = Path(in_dir), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    hh_book = Path(hh_path) if hh_path else None

    panels = [
        ("A", "eGFP + HEPES"),
        ("B", "eGFP + TRIS"),
        ("C", "CSP + HEPES"),
        ("D", "CSP + TRIS"),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(BASE_W, BASE_H), sharex=True, sharey=True)

    for ax, (letter, pretty) in zip(axes.flatten(), panels):
        sheet = SHEET_MAP[pretty]
        # experimental RNA from unified Excel (rna2)
        t_exp, y_exp = load_exp_rna(Path(all_data_path), sheet)
        # optional HH overlay (legacy single-sheet or multi-sheet)
        hh_vals = load_hh_series(hh_book, sheet) if hh_book else None
        # UKF CSV path
        ukf_csv = in_dir / CSV_DEFAULTS[pretty]
        plot_panel(ax, ukf_csv, t_exp, y_exp, hh_vals, letter,
                   show_uncertainty=show_uncertainty,
                   unc_scale=(unc if isinstance(unc, (int, float)) else None),
                   smooth=smooth, window=window,
                   plot_model_at_exp=plot_model_at_exp)
        # title text after panel letter
        ax.text(0.09, PANEL_Y, pretty, transform=ax.transAxes, fontsize=LABFS, va='bottom', ha='left')

    # legend (unique labels)
    handles, labels = [], []
    for ax in axes.flatten():
        h, l = ax.get_legend_handles_labels()
        for hh, ll in zip(h, l):
            if ll not in labels:
                handles.append(hh); labels.append(ll)
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, LEG_Y),
               ncol=min(4, len(labels)), frameon=False, fontsize=LABFS)

    fig.subplots_adjust(left=0.10, right=0.995, bottom=0.12, top=TOP, wspace=WSPACE, hspace=HSPACE)
    fig.savefig(Path(out_dir) / "quad_rna_profiles.pdf", bbox_inches='tight')
    fig.savefig(Path(out_dir) / "quad_rna_profiles.png", dpi=600, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    hh_auto = resolve_hh_path()
    quad_plot(
        in_dir=UKF_DIR,
        all_data_path=ALL_DATA_XLSX,
        hh_path=hh_auto,         
        out_dir=OUT_DIR,
        show_uncertainty=False,   # set True to show ±2σ bands
        unc=1.0,                  
        smooth=True,
        window=20,
        plot_model_at_exp=False   # True = draw UKF at experimental timestamps
    )
