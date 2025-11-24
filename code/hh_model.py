#!/usr/bin/env python3
"""
hh_model.py
===========
Henderson–Hasselbalch (HH) utilities for multi-buffer mixtures used in IVT / bioprocess media.

Author: Mahdi Ahmed, Shady Hamed, 2025
License: MIT
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable, Tuple, Dict, List

import numpy as np
import pandas as pd
from scipy.optimize import brentq
import matplotlib

import matplotlib.pyplot as plt


LOG10 = np.log(10.0)
KW = 1.0e-14  # water ion product at 25°C (approx)

@dataclass(frozen=True)
class Buffer:
    """Monoprotic buffer specification."""
    name: str
    pKa: float        # acid dissociation constant (pKa)
    C: float          # total analytical concentration (mol/L)

    @property
    def Ka(self) -> float:
        return 10.0 ** (-self.pKa)


# ---------------------------- Core HH formulas ----------------------------

def species_distribution(buffer: Buffer, pH: float) -> Tuple[float, float]:
    """
    Return (HA, A-) concentrations at given pH for a single monoprotic buffer.
    Uses fractions HA = C * H / (Ka + H), A- = C * Ka / (Ka + H).
    """
    H = 10.0 ** (-pH)
    Ka = buffer.Ka
    denom = (Ka + H)
    HA = buffer.C * H / denom
    A_ = buffer.C * Ka / denom
    return HA, A_


def mixture_distribution(buffers: Iterable[Buffer], pH: float) -> pd.DataFrame:
    """
    Compute species distributions for all buffers at target pH; returns a DataFrame with
    columns: name, pKa, C, H, Ka, HA, A_minus, alpha_acid, alpha_base.
    """
    H = 10.0 ** (-pH)
    rows = []
    for b in buffers:
        Ka = b.Ka
        HA, A_ = species_distribution(b, pH)
        rows.append({
            "name": b.name, "pKa": b.pKa, "C": b.C, "H": H, "Ka": Ka,
            "HA": HA, "A_minus": A_,
            "alpha_acid": HA / max(b.C, 1e-30),
            "alpha_base": A_ / max(b.C, 1e-30),
        })
    df = pd.DataFrame(rows).set_index("name")
    df.attrs["pH"] = pH
    return df


def buffer_capacity(buffers: Iterable[Buffer], pH: float, include_water: bool = True) -> float:
    """
    Total buffer capacity β = dB/dpH (base added per pH unit), for a mixture of monoprotic buffers.
    For each acid: β_i = 2.303 * C_i * Ka*H/(Ka + H)^2 . Optionally add water: β_w = 2.303*(H + Kw/H).
    """
    H = 10.0 ** (-pH)
    beta = 0.0
    for b in buffers:
        Ka = b.Ka
        beta += (2.303) * b.C * (Ka * H) / ((Ka + H) ** 2)
    if include_water:
        beta += (2.303) * (H + (KW / H))
    return float(beta)


# ---------------------------- Titration utilities ----------------------------

def grid_summary(buffers: Iterable[Buffer], pH_min: float, pH_max: float, points: int = 1001) -> pd.DataFrame:
    """
    Compute summary across a pH grid: total A-, total HA, capacity.
    Returns a tidy DataFrame with index pH and columns:
        total_A_minus, total_HA, beta, and per-component A_minus_<name> if unique names.
    """
    pH = np.linspace(pH_min, pH_max, points)
    H = 10.0 ** (-pH)
    buffers = list(buffers)
    out = {
        "pH": pH,
        "total_A_minus": np.zeros_like(pH),
        "total_HA": np.zeros_like(pH),
        "beta": np.zeros_like(pH),
    }
    comp_cols: Dict[str, np.ndarray] = {}
    for b in buffers:
        Ka = b.Ka
        denom = (Ka + H)
        A_ = b.C * Ka / denom
        HA = b.C * H / denom
        out["total_A_minus"] += A_
        out["total_HA"] += HA
        comp_cols[f"Aminus_{b.name}"] = A_
    out["beta"] = np.array([buffer_capacity(buffers, x, include_water=True) for x in pH])
    df = pd.DataFrame(out).set_index("pH")
    for k, v in comp_cols.items():
        df[k] = v
    return df


def plot_systems(systems: List[Tuple[str, List[Buffer]]],
                 pH_range: Tuple[float, float] = (3.0, 9.0),
                 points: int = 1000,
                 outdir: str | Path = "figures/hh",
                 show: bool = False) -> List[Path]:
    """
    Plot total A- vs pH (and capacity on secondary axis) for a list of systems.
    Each 'system' is (title, [Buffer(...), ...]).
    Saves PNG files under outdir and returns their paths.
    """
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    paths: List[Path] = []
    pH_vals = np.linspace(pH_range[0], pH_range[1], points)
    cols = 2
    rows = int(np.ceil(len(systems) / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(12, 4 * rows), squeeze=False)

    for idx, (title, bufs) in enumerate(systems):
        r, c = divmod(idx, cols)
        ax = axes[r][c]
        # total A- and optionally components
        df = grid_summary(bufs, *pH_range, points=points)
        ax.plot(df.index, df["total_A_minus"], lw=2, label="Σ A⁻")
        # plot components faintly
        for b in bufs:
            ax.plot(df.index, df[f"Aminus_{b.name}"], lw=1, alpha=0.6, label=f"A⁻ {b.name}")
        ax.set_title(title, loc="left", fontweight="bold", pad=6)
        ax.set_xlabel("pH"); ax.set_ylabel("Total [A⁻] (mol/L)")
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize="small", ncol=2, loc="upper left")

        # Secondary axis for buffer capacity
        ax2 = ax.twinx()
        ax2.plot(df.index, df["beta"], lw=1.5, ls="--", alpha=0.75)
        ax2.set_ylabel("β (dB/dpH)")

    for j in range(len(systems), rows * cols):
        r, c = divmod(j, cols)
        fig.delaxes(axes[r][c])

    fig.tight_layout()
    out = outdir / "hh_systems_summary.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    else:
        plt.close(fig)
    paths.append(out)
    return paths


def calculate_initial_concentrations(
    buffers: Iterable[Buffer], target_pH: float
) -> pd.DataFrame:
    """
    Convenience wrapper matching the intent of the notebook's
    'calculate_initial_concentrations(...)' but generalized.

    Parameters
    ----------
    buffers : iterable of Buffer
    target_pH : float

    Returns
    -------
    DataFrame with per-buffer HA, A-, alpha fractions, and totals stored in attrs.
    """
    df = mixture_distribution(buffers, target_pH)
    df.attrs["total_A_minus"] = float(df["A_minus"].sum())
    df.attrs["total_HA"] = float(df["HA"].sum())
    df.attrs["beta"] = buffer_capacity(buffers, target_pH, include_water=True)
    return df


def default_systems() -> List[Tuple[str, List[Buffer]]]:
    """Replicate the co-author's grid of systems (HEPES vs TRIS at various levels)."""
    def B(name, pKa, C): return Buffer(name=name, pKa=pKa, C=C)
    # Fixed components
    acetate = ("Acetate", 4.73, 0.042)  # or 0.082 when TRIS acetate used
    ntp     = ("NTP", 6.95, 0.01)
    phosphate = ("Pi", 7.20, 0.06)
    hepes = ("HEPES", 7.30)
    tris  = ("TRIS", 7.80)

    systems: List[Tuple[str, List[Buffer]]] = [
        ("A) HEPES 20 mM",
         [B(*acetate), B(*ntp), B(*phosphate), B(hepes[0], hepes[1], 0.02)]),
        ("B) TRIS acetate 20 mM",
            [B("Acetate", 4.73, 0.082), B(*ntp), B(*phosphate), B(tris[0], tris[1], 0.02)]),
        ("C) HEPES 40 mM",
         [B(*acetate), B(*ntp), B(*phosphate), B(hepes[0], hepes[1], 0.04)]),
        ("D) TRIS acetate 40 mM",
            [B("Acetate", 4.73, 0.082), B(*ntp), B(*phosphate), B(tris[0], tris[1], 0.04)]),
        ("E) HEPES 80 mM",
         [B(*acetate), B(*ntp), B(*phosphate), B(hepes[0], hepes[1], 0.08)]),
        ("F) TRIS acetate 80 mM",
            [B("Acetate", 4.73, 0.082), B(*ntp), B(*phosphate), B(tris[0], tris[1], 0.08)]),
        ("G) HEPES 160 mM",
         [B(*acetate), B(*ntp), B(*phosphate), B(hepes[0], hepes[1], 0.16)]),
        ("H) TRIS acetate 160 mM",
            [B("Acetate", 4.73, 0.082), B(*ntp), B(*phosphate), B(tris[0], tris[1], 0.16)]),
    ]
    return systems


def demo(outdir: str | Path = "figures/hh", show: bool = False) -> List[Path]:
    """Run the default plots and save a CSV snapshot of one system at pH 7.0."""
    # 1) Plots for all systems
    paths = plot_systems(default_systems(), outdir=outdir, show=show)

    title, bufs = default_systems()[2]
    df = calculate_initial_concentrations(bufs, target_pH=7.0)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "distribution_HEPES_40mM_pH7.0.csv"
    df.to_csv(csv_path)
    paths.append(csv_path)
    return paths


def calculate_initial_concentrations_legacy(
    hepes_pka: float, hepes_conc: float,
    acetate_pka: float, acetate_conc: float,
    pi_pka: float, pi_conc: float,
    ntp_pka: float, ntp_conc: float,
    target_pH: float
) -> pd.DataFrame:
    """
    Backwards-compatible helper to mirror a prior notebook signature that passed
    four (pKa, C) pairs explicitly (HEPES, acetate, Pi, NTP) plus target_pH.

    Returns a DataFrame like mixture_distribution(...).
    """
    bufs = [
        Buffer("HEPES", hepes_pka, hepes_conc),
        Buffer("Acetate", acetate_pka, acetate_conc),
        Buffer("Pi", pi_pka, pi_conc),
        Buffer("NTP", ntp_pka, ntp_conc),
    ]
    return calculate_initial_concentrations(bufs, target_pH)


if __name__ == "__main__":
    outs = demo(outdir=Path("figures") / "hh_demo", show=False)
    for p in outs:
        print(f"[OK] wrote {p}")
