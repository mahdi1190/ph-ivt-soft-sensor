#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot_species_grid.py
====================

Grid figure for selected IVT species across four datasets (eGFP/CSP × HEPES/TRIS).

Inputs (default locations):
  reports/softsensor_run/
    - ivt_ukf_results_egfp_HEPES.csv
    - ivt_ukf_results_egfp_TRIS.csv
    - ivt_ukf_results_csp_HEPES.csv
    - ivt_ukf_results_csp_TRIS.csv

Outputs:
  figures/species/superplot.pdf
  figures/species/superplot.png

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
from matplotlib.ticker import ScalarFormatter

# ---------- unified layout (7.3 × 8.0, compact, all fonts = 10) ----------
FIGSIZE         = (7.3, 8.0)
TOP_FIG         = 0.80
BOTTOM_FIG      = 0.10
LEFT_FIG        = 0.11
RIGHT_FIG       = 0.995
WSPACE          = 0.34
HSPACE          = 0.48

# Legend placement
LEGEND_Y1       = 0.90

# Typography (all 10)
AX_LABEL_SIZE   = 10
TICK_LABEL_SIZE = 10
PANEL_LABEL_SIZE= 10
LEGEND_FS       = 10

# Pads
XLABEL_PAD      = 2
YLABEL_PAD      = 2
XTICK_PAD       = 2
YTICK_PAD       = 2

# Lines
LINEWIDTH       = 1.5

# Panel headers
PANEL_TITLE_Y   = 1.030
TITLE_X_OFFSET  = 0.10

CSV_DEFAULTS = {
    "eGFP + HEPES": "ivt_ukf_results_egfp_HEPES.csv",
    "eGFP + TRIS":  "ivt_ukf_results_egfp_TRIS.csv",
    "CSP + HEPES":  "ivt_ukf_results_csp_HEPES.csv",
    "CSP + TRIS":   "ivt_ukf_results_csp_TRIS.csv",
}

# Turn off scientific notation and offsets globally for axes
plt.rcParams['axes.formatter.useoffset'] = False
plt.rcParams['axes.formatter.use_locale'] = False

def pick_column(df: pd.DataFrame, include_keys=(), exclude_keys=()) -> str | None:
    """Token-based fuzzy matcher: all include tokens present, no exclude tokens."""
    include = [k.lower() for k in include_keys]
    exclude = [k.lower() for k in exclude_keys]
    for col in df.columns:
        tokens = set(t for t in re.split(r'[^a-z0-9]+', col.lower()) if t)
        if all(k in tokens for k in include) and not any(x in tokens for x in exclude):
            return col
    return None

def numeric_series(df: pd.DataFrame, col: str) -> pd.Series:
    return pd.to_numeric(df[col], errors="coerce")

def finite_mask(*arrs) -> np.ndarray:
    m = np.ones_like(np.asarray(arrs[0], dtype=float), dtype=bool)
    for a in arrs:
        m &= np.isfinite(np.asarray(a, dtype=float))
    return m

def units_and_factor(species: str) -> tuple[str, float]:
    """
    CSV values are in mol/L; convert to nicer display units.
    """
    s = species.lower()
    if s in {"ppi", "pi", "mgatp", "mgctp", "mggtp", "mgutp"}:
        return "mM", 1e6   # show as mM
    return "µM", 1e6       # fallback as µM

def apply_plain_ticks(ax: plt.Axes):
    """Force plain numeric ticks (no 1e-3 multipliers, no offsets)."""
    fmt = ScalarFormatter(useOffset=False, useMathText=False)
    fmt.set_scientific(False)
    ax.yaxis.set_major_formatter(fmt)
    # (optional) also for x if needed:
    # ax.xaxis.set_major_formatter(fmt)

def plot_species_grid(
    dataset_paths: dict[str, str | Path],
    species_list: list[str],
    smooth: bool = False,
    window: int = 5,
    cols: int = 2,
    figsize: tuple = FIGSIZE,
    out_path: str | Path = "figures/species/superplot.pdf",
    also_png: bool = True,
):
    # styles per dataset
    palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
    styles = {
        label: dict(color=palette[i % len(palette)], linestyle='-', lw=LINEWIDTH)
        for i, label in enumerate(dataset_paths)
    }

    n = len(species_list)
    rows = int(np.ceil(n / cols))
    fig, axes_arr = plt.subplots(rows, cols, figsize=figsize)
    axes = np.ravel(axes_arr)

    for idx, species in enumerate(species_list):
        ax = axes[idx]
        unit_label, factor = units_and_factor(species)

        for label, path in dataset_paths.items():
            df = pd.read_csv(path)
            # time
            tcol = pick_column(df, include_keys=('time','min')) or 'time_min'
            t = numeric_series(df, tcol)

            # values
            ycol = pick_column(df, include_keys=(species.lower(),), exclude_keys=('minus','plus')) or species
            if ycol not in df.columns:
                raise KeyError(f"{species}: could not find a column like '{species}' in {path}")
            y = numeric_series(df, ycol) * factor
            if smooth:
                y = y.rolling(window=window, center=True, min_periods=1).mean()
            y = y.clip(lower=0)

            m = finite_mask(t, y)
            if m.any():
                ax.plot(t[m], y[m], label=label, **styles[label])

        # y-limits and tick formatting
        if species.lower() != 'ph_pred':
            ax.set_ylim(bottom=0)
        apply_plain_ticks(ax)  # <- no scientific 1e−3, no offsets

        # Panel header
        letter = chr(65 + idx)
        ax.text(0.0, PANEL_TITLE_Y, f'({letter})',
                transform=ax.transAxes, ha='left', va='bottom',
                fontsize=PANEL_LABEL_SIZE, fontweight='bold', clip_on=False)
        ax.text(TITLE_X_OFFSET, PANEL_TITLE_Y, species,
                transform=ax.transAxes, ha='left', va='bottom',
                fontsize=PANEL_LABEL_SIZE, fontweight='normal', clip_on=False)

        # Axes cosmetics (unified)
        ax.set_xlabel('Time (minutes)', fontsize=AX_LABEL_SIZE, labelpad=XLABEL_PAD)
        ax.set_ylabel(f'Concentration ({unit_label})', fontsize=AX_LABEL_SIZE, labelpad=YLABEL_PAD)
        ax.set_xlim(0, 120)
        ax.set_xticks(np.arange(0, 121, 20))
        ax.tick_params(axis='x', labelsize=TICK_LABEL_SIZE, pad=XTICK_PAD)
        ax.tick_params(axis='y', labelsize=TICK_LABEL_SIZE, pad=YTICK_PAD)
        ax.grid(False)

    # Remove extra axes if grid > n
    for ax in axes[n:]:
        ax.remove()

    # ---------- Dataset legend (2 columns → 2 rows for 4 items) ----------
    ds_handles = [Line2D([0],[0], color=styles[l]['color'],
                         lw=styles[l]['lw'], linestyle=styles[l]['linestyle'])
                  for l in dataset_paths]
    ds_labels  = list(dataset_paths.keys())
    fig.legend(handles=ds_handles, labels=ds_labels,
               loc='upper center', ncol=2, frameon=False, fontsize=LEGEND_FS,
               bbox_to_anchor=(0.5, LEGEND_Y1),
               borderaxespad=0.0, handletextpad=0.6,
               columnspacing=1.2, labelspacing=0.3)

    # Final layout & export
    fig.subplots_adjust(left=LEFT_FIG, right=RIGHT_FIG,
                        bottom=BOTTOM_FIG, top=TOP_FIG,
                        wspace=WSPACE, hspace=HSPACE)

    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=600, bbox_inches='tight')
    if also_png:
        fig.savefig(out_path.with_suffix(".png"), dpi=600, bbox_inches='tight')
    plt.close(fig)

# --------- CLI ---------
def main():
    ap = argparse.ArgumentParser(description="Grid plot of IVT species across four runs (no CI bands, plain ticks).")
    ap.add_argument("--in-dir", default="reports/softsensor_run",
                    help="Directory containing the four UKF CSVs.")
    ap.add_argument("--out", default="figures/species/superplot.pdf",
                    help="Output PDF path (PNG written alongside).")
    ap.add_argument("--species", nargs="+",
                    default=["PPi","Pi","MgATP","MgCTP","MgGTP","MgUTP"],
                    help="Species names as they appear (prefix) in CSV headers.")
    ap.add_argument("--smooth", action="store_true", help="Enable rolling mean smoothing.")
    ap.add_argument("--window", type=int, default=12, help="Smoothing window size.")
    ap.add_argument("--cols", type=int, default=2, help="Number of columns in the grid.")
    args = ap.parse_args()

    in_dir = Path(args.in_dir)
    dataset_paths = {
        "eGFP + HEPES": in_dir / CSV_DEFAULTS["eGFP + HEPES"],
        "eGFP + TRIS":  in_dir / CSV_DEFAULTS["eGFP + TRIS"],
        "CSP + HEPES":  in_dir / CSV_DEFAULTS["CSP + HEPES"],
        "CSP + TRIS":   in_dir / CSV_DEFAULTS["CSP + TRIS"],
    }

    plot_species_grid(
        dataset_paths=dataset_paths,
        species_list=args.species,
        smooth=args.smooth,
        window=int(args.window),
        cols=int(args.cols),
        figsize=FIGSIZE,
        out_path=args.out,
        also_png=True,
    )
    print(f"[OK] wrote {args.out} and {Path(args.out).with_suffix('.png')}")

if __name__ == "__main__":
    main()
