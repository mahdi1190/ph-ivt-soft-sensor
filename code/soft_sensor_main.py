#!/usr/bin/env python3
"""
soft_sensor_main.py
===================

Integrated IVT DAE model (36 states) + Unscented Kalman Filter (UKF) soft sensor for pH,
packaged so you can either import a single function or run this file directly.

This module exposes:
- default_kin_params()         → kinetic parameter dict
- make_initials(pH_target)     → consistent 36-state initial vector
- run_ivt(...)                 → integrates the DAE (Assimulo/IDA)
- run_ukf(...)                 → UKF over the DAE with robust PD covariance handling
- load_and_sample_ph(path, ...)→ helper to load CSV/Excel pH data and sample uniformly
- run_soft_sensor(...)         → ONE CALL: load → UKF → CSV + figures

Notes
-----
• Requires Assimulo (with SUNDIALS). Recommended:
    conda install -c conda-forge assimulo sundials
• Figures are saved to files (non-interactive backend used if needed).
• Constants match model (your K_MgNTP (K1) and K_DNAT7 (K12) preserved).

Author: Your Name (2025)
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Any, Literal, Tuple, List

import numpy as np
import pandas as pd
from scipy.optimize import root_scalar

# --- plotting (headless-safe) ---
import matplotlib
if not matplotlib.get_backend().lower().startswith("qt"):
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Attempt Assimulo import (friendly error if missing)
try:
    from assimulo.problem import Implicit_Problem
    from assimulo.solvers import IDA
except Exception as exc:  # pragma: no cover
    raise RuntimeError(
        "Assimulo with SUNDIALS is required. Install via conda:\n"
        "  conda install -c conda-forge assimulo sundials\n"
        "or via pip (ensure system SUNDIALS is available)."
    ) from exc

# UKF (FilterPy)
from filterpy.kalman import UnscentedKalmanFilter, MerweScaledSigmaPoints

# -------------------------------------------------------------------------
# 1) CONSTANTS
# -------------------------------------------------------------------------
KW      = 1.0e-14              # water autodissociation
K_HEPES = 10.0 ** -8.3         # HEPES Ka
K_ACET  = 10.0 ** -4.76        # acetate Ka
K_PI    = 10.0 ** -6.15        # H2PO4- ⇌ H+ + HPO4^2-
K_NTP   = 10.0 ** -6.6         # simplified NTP buffering Ka

# Mg–NTP association “K” (effective) — your values
K_MgNTP: Dict[str, float] = {"ATP": 0.0005, "GTP": 0.0005, "CTP": 0.0005, "UTP": 0.0005}

# Your values
K_DNAT7    = 110.0             # µM-scale DNA·T7 KD
KSP_MG2PPi = 1.4e-6            # [PPi]eq at 1 M Mg²⁺
UNIT_ACTIVITY = 16.67e-6       # mol L⁻¹ min⁻¹
SPEC_ACTIVITY = 2.0e8          # Units per mol RNAP

# Exact nucleotide counts for the transcript template
NA, NC, NG, NU = 266, 279, 243, 142
NUC_COUNTS = {"ATP": NA, "CTP": NC, "GTP": NG, "UTP": NU}
N_TOTAL = NA + NC + NG + NU

# -------------------------------------------------------------------------
# 2) KINETIC PARAMETERS
# -------------------------------------------------------------------------
def default_kin_params() -> Dict[str, float]:
    """Default kinetic/side-process parameters."""
    return {
        "k_app": 1.0,          # min⁻¹
        "k1": 1.0e4,           # Mg²⁺ inhibition coefficient
        "k2": 1.0e-3,          # MgNTP inhibition coefficient
        "k_PPase": 3120.0,     # PPi → 2 Pi hydrolysis
        "k_prec": 500.0,       # Mg₂PPi precipitation
        "k_d": 0.0,            # RNAP decay (disabled if zero)
        "k_ds": 0.0,           # dsRNA mis-anneal
        "init_seq_len": 8,     # nucleotides
        "init_rate_nt_s": 0.5, # nt/s
        "elong_rate_nt_s": 200.0,  # nt/s
    }

# -------------------------------------------------------------------------
# 3) STATE INDEXING (36 states) + ALGEBRAIC FLAGS
# -------------------------------------------------------------------------
RNA, PPi, Pi = 0, 1, 2
ATPt, CTPt, GTPt, UTPt = 3, 4, 5, 6
Htot, Mgtot = 7, 8
DNAt, T7tot = 9, 10
PPis, dsRNA = 11, 12
HEPtot, ACtot, SPDtot = 13, 14, 15
Hfree, Mgfree = 16, 17
MgATP, MgCTP, MgGTP, MgUTP = 18, 19, 20, 21
HEP, HHEP = 22, 23
Ac, HAc = 24, 25
H2Pi, HPi = 26, 27
ATP_base, HATP = 28, 29
CTP_base, HCTP = 30, 31
GTP_base, HGTP = 32, 33
UTP_base, HUTP = 34, 35

N_STATE = 36
# differential (1) vs algebraic (0)
ALGVAR = np.array([1]*16 + [0]*(N_STATE-16), dtype=int)

# For nicer outputs/plots
STATE_NAMES: List[str] = [
    "RNA", "PPi", "Pi",
    "ATP_tot", "CTP_tot", "GTP_tot", "UTP_tot",
    "H_tot", "Mg_tot",
    "DNA_tot", "T7_tot",
    "PPi_solid", "dsRNA",
    "HEPES_tot", "Acetate_tot", "Spd_tot",
    "H_free", "Mg_free",
    "MgATP", "MgCTP", "MgGTP", "MgUTP",
    "HEP_base", "HEP_acid",
    "Ac_base", "Ac_acid",
    "H2Pi", "HPi",
    "ATP_base", "HATP",
    "CTP_base", "HCTP",
    "GTP_base", "HGTP",
    "UTP_base", "HUTP",
]

# -------------------------------------------------------------------------
# 4) INITIAL CONDITIONS
# -------------------------------------------------------------------------
def default_initials(pH_target: float = 6.67) -> np.ndarray:
    """Base initial state with explicit buffer/NTP splits (not charge-consistent yet)."""
    y = np.zeros(N_STATE)

    # main reactant pools
    y[ATPt:UTPt+1] = 0.010     # M each free nucleotide
    y[Mgtot] = 0.042           # M
    y[DNAt]  = 50.0            # µM
    y[T7tot] = 246.1           # Units/µL
    y[HEPtot] = 0.040          # mM
    y[ACtot]  = 0.042
    y[SPDtot] = 1e-12          # M
    y[Pi]     = 1e-12
    y[PPis]   = 1e-12
    y[dsRNA]  = 1e-12

    # target free protons
    H = 10**(-pH_target)
    y[Hfree] = H

    # buffer splits
    y[HEP]  = y[HEPtot] * K_HEPES / (K_HEPES + H)
    y[HHEP] = y[HEPtot] - y[HEP]

    y[Ac]   = y[ACtot] * K_ACET / (K_ACET + H)
    y[HAc]  = y[ACtot] - y[Ac]

    # phosphate split (H2Pi acid, HPi base)
    y[H2Pi] = y[Pi] / (1 + K_PI / H)
    y[HPi]  = y[Pi] - y[H2Pi]

    # NTP buffering splits
    for base, idx_base, idx_acid, idx_tot in [
        ("ATP", ATP_base, HATP, ATPt),
        ("CTP", CTP_base, HCTP, CTPt),
        ("GTP", GTP_base, HGTP, GTPt),
        ("UTP", UTP_base, HUTP, UTPt),
    ]:
        y[idx_base] = y[idx_tot] * K_NTP / (K_NTP + H)
        y[idx_acid] = y[idx_tot] - y[idx_base]

    # approx proton pool
    y[Htot] = H + KW / H + y[HHEP] + y[HAc] + y[H2Pi] + y[HATP] + y[HCTP] + y[HGTP] + y[HUTP]

    # Mg guesses
    y[Mgfree] = 0.0
    y[MgATP:MgUTP+1] = 0.0
    return y


def make_initials(pH_target: float = 6.67) -> np.ndarray:
    """Adjust H_tot such that free [H+] = 10^(-pH_target) after algebraic equilibria."""
    y = default_initials(pH_target)
    H_target = 10 ** (-pH_target)
    y[Hfree] = H_target

    def charge_residual(H_tot: float) -> float:
        hep_acid = y[HEPtot] * H_target / (K_HEPES + H_target)
        ac_acid  = y[ACtot] * H_target / (K_ACET + H_target)
        pi_acid  = y[Pi] / (1 + K_PI / H_target)  # H2Pi
        atp_acid = y[ATPt] * H_target / (K_NTP + H_target)
        ctp_acid = y[CTPt] * H_target / (K_NTP + H_target)
        gtp_acid = y[GTPt] * H_target / (K_NTP + H_target)
        utp_acid = y[UTPt] * H_target / (K_NTP + H_target)
        return (
            H_tot - (H_target + KW / H_target) - hep_acid - ac_acid - pi_acid
            - atp_acid - ctp_acid - gtp_acid - utp_acid
        )

    hi = max(y[Htot], 0.1)
    sol = root_scalar(charge_residual, bracket=[1e-20, 10 * hi], method="bisect")
    if not sol.converged:  # pragma: no cover
        raise RuntimeError("Failed to solve for consistent H_tot in make_initials().")
    y[Htot] = sol.root
    return y

# -------------------------------------------------------------------------
# 5) HELPERS
# -------------------------------------------------------------------------
def geom_mean_lim(mgntp: Dict[str, float]) -> float:
    """Geometric-mean limiter over MgNTP species, weighted by nucleotide counts."""
    prod = 1.0
    for base, ν in NUC_COUNTS.items():
        prod *= mgntp[base] / ν
    return prod ** (1 / len(NUC_COUNTS))

# -------------------------------------------------------------------------
# 6) DAE RESIDUAL
# -------------------------------------------------------------------------
def residual(t: float, y: np.ndarray, yd: np.ndarray, p: Dict[str, float]) -> np.ndarray:
    """DAE residual g(t, y, ydot) = 0 for the 36-state system."""
    res = np.zeros_like(y)
    H, Mg = y[Hfree], y[Mgfree]

    # (A) explicit buffer algebra (HEPES, acetate, phosphate)
    HHEP_calc = H * y[HEP] / K_HEPES
    HAc_calc  = H * y[Ac]  / K_ACET
    H2Pi_calc = H * y[HPi] / K_PI

    res[HEP]  = y[HEP]  + y[HHEP] - y[HEPtot]
    res[HHEP] = y[HHEP] - HHEP_calc
    res[Ac]   = y[Ac]   + y[HAc]  - y[ACtot]
    res[HAc]  = y[HAc]  - HAc_calc

    res[HPi]  = y[HPi] + y[H2Pi] - y[Pi]
    res[H2Pi] = y[H2Pi] - H2Pi_calc

    # (B) NTP buffering algebra
    def ntp_pair(base_idx: int, acid_idx: int, tot_idx: int) -> None:
        Hacid_calc = H * y[base_idx] / K_NTP
        res[base_idx] = y[base_idx] + y[acid_idx] - y[tot_idx]
        res[acid_idx] = y[acid_idx] - Hacid_calc

    ntp_pair(ATP_base, HATP, ATPt)
    ntp_pair(CTP_base, HCTP, CTPt)
    ntp_pair(GTP_base, HGTP, GTPt)
    ntp_pair(UTP_base, HUTP, UTPt)

    # Proton balance
    res[Hfree] = (
        y[Htot] - (H + KW / H) - y[HHEP] - y[HAc] - y[H2Pi] - y[HATP] - y[HCTP] - y[HGTP] - y[HUTP]
    )

    # (C) Mg–NTP complexes
    mgntp = {"ATP": y[MgATP], "CTP": y[MgCTP], "GTP": y[MgGTP], "UTP": y[MgUTP]}
    ntp_free = {
        "ATP": y[ATPt] - mgntp["ATP"],
        "CTP": y[CTPt] - mgntp["CTP"],
        "GTP": y[GTPt] - mgntp["GTP"],
        "UTP": y[UTPt] - mgntp["UTP"],
    }
    for base, K in K_MgNTP.items():
        idx = {"ATP": MgATP, "CTP": MgCTP, "GTP": MgGTP, "UTP": MgUTP}[base]
        res[idx] = mgntp[base] - (Mg * ntp_free[base] / K)

    res[Mgfree] = y[Mgtot] - (Mg + sum(mgntp.values()) + 2.0 * y[PPis])

    # (D) transcription & side-rates
    lim_sub = 0.0 if min(mgntp.values()) < 1e-12 else geom_mean_lim(mgntp)
    DNAT7 = (y[DNAt] * y[T7tot]) / K_DNAT7
    V_tr_raw = (
        p["k_app"] * Mg * DNAT7 * lim_sub /
        (1 + p["k1"] * Mg + p["k2"] * sum(mgntp.values()))
    )

    # capacity terms (kept for context; no ceiling applied in this variant)
    L = N_TOTAL
    τ_init = p["init_seq_len"] / (p["init_rate_nt_s"] * 60)
    _throughput = min(1 / τ_init, p["elong_rate_nt_s"] * 60 / L)
    V_tr = V_tr_raw

    V_prec  = p["k_prec"] * max(0.0, y[PPi] - KSP_MG2PPi / (Mg**2 + 1e-30))
    V_ppase = p["k_PPase"] * y[PPi] / (1 + y[PPi])
    V_ds    = p["k_ds"]   * y[RNA] * Mg

    νA, νC, νG, νU = NUC_COUNTS.values()

    # (E) differential balances
    res[RNA]   = yd[RNA]  - (V_tr - V_ds)
    res[PPi]   = yd[PPi]  - (N_TOTAL * V_tr - V_ppase - V_prec)
    res[Pi]    = yd[Pi]   - (2.0 * V_ppase)

    res[ATPt]  = yd[ATPt] + νA * V_tr
    res[CTPt]  = yd[CTPt] + νC * V_tr
    res[GTPt]  = yd[GTPt] + νG * V_tr
    res[UTPt]  = yd[UTPt] + νU * V_tr

    res[Htot]  = yd[Htot] - (N_TOTAL * V_tr)
    res[Mgtot] = yd[Mgtot] - 2.0 * V_prec

    res[DNAt]  = yd[DNAt]         # constant
    res[T7tot] = yd[T7tot]        # constant

    res[PPis]  = yd[PPis] - V_prec
    res[dsRNA] = yd[dsRNA] - V_ds

    # totals held constant
    res[HEPtot] = yd[HEPtot]
    res[ACtot]  = yd[ACtot]
    res[SPDtot] = yd[SPDtot]

    return res

# -------------------------------------------------------------------------
# 7) IVT INTEGRATOR
# -------------------------------------------------------------------------
def run_ivt(
    t_final: float = 1.0,
    ncp: int = 1,
    initials: np.ndarray | None = None,
    kinpars: Dict[str, float] | None = None,
    atol: float = 1e-12,
    rtol: float = 1e-3,
):
    """Integrate the DAE with Assimulo/IDA and return (t, y, yd)."""
    y0 = make_initials(pH_target=6.67) if initials is None else initials.copy()
    p = kinpars if kinpars is not None else default_kin_params()
    yd0 = np.zeros_like(y0)
    problem = Implicit_Problem(lambda t, y, yd: residual(t, y, yd, p), y0, yd0)
    problem.algvar = ALGVAR
    solver = IDA(problem)
    solver.atol = atol
    solver.rtol = rtol
    solver.make_consistent('IDA_YA_YDP_INIT')
    return solver.simulate(t_final, ncp=ncp)

# -------------------------------------------------------------------------
# 8) UKF helpers
# -------------------------------------------------------------------------
def _propagate_state(
    y0: np.ndarray,
    dt: float,
    params: Dict[str, float],
    atol: float = 1e-12,
    rtol: float = 1e-3,
) -> np.ndarray:
    yd0 = np.zeros_like(y0)
    prob = Implicit_Problem(lambda t, y, yd: residual(t, y, yd, params), y0, yd0)
    prob.algvar = ALGVAR
    solver = IDA(prob)
    solver.atol = atol
    solver.rtol = rtol
    solver.make_consistent('IDA_YA_YDP_INIT')
    _, y_arr, _ = solver.simulate(dt, ncp=1)
    return y_arr[-1].copy()

def _fx(y: np.ndarray, dt: float, params: Dict[str, float]) -> np.ndarray:
    return _propagate_state(y, dt, params)

def _hx(y: np.ndarray) -> np.ndarray:
    h_conc = max(y[Hfree], 1e-20)
    ph_pred = -np.log10(h_conc)
    return np.array([ph_pred])

def _make_pos_def(A: np.ndarray, epsilon: float = 1e-12) -> np.ndarray:
    """Force symmetric PD with eigenvalue floor epsilon."""
    A = (A + A.T) / 2
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.clip(eigvals, epsilon, None)
    return eigvecs @ np.diag(eigvals) @ eigvecs.T

# -------------------------------------------------------------------------
# 9) UKF DRIVER
# -------------------------------------------------------------------------
def run_ukf(
    t_meas: np.ndarray,            # shape (N,)
    ph_meas: np.ndarray,           # shape (N,)
    y0: np.ndarray,
    params: Dict[str, float],
    r_std: float = 1e-4,
    alpha: float = 1e-3,
    beta: float = 2.0,
    kappa: float = 0.0,
    noise_mode: Literal["noisy","smooth"] = "noisy",
    verbose: bool = False,
) -> Dict[str, Any]:
    """Unscented Kalman Filter over the full 36-state DAE."""
    if t_meas.shape != ph_meas.shape:
        raise ValueError("t_meas and ph_meas must have the same length.")

    dt_seq = np.diff(t_meas, prepend=t_meas[0])

    # Sigma points + UKF
    sigmas = MerweScaledSigmaPoints(n=N_STATE, alpha=alpha, beta=beta, kappa=kappa)
    ukf = UnscentedKalmanFilter(
        dim_x=N_STATE, dim_z=1,
        fx=lambda x, dt: _fx(x, dt, params), hx=_hx, dt=1.0, points=sigmas,
    )

    # Indices most coupled to pH dynamics
    coupled_idx = {Hfree, Htot, Mgfree, HHEP, HAc, MgATP, MgCTP, MgGTP, MgUTP}

    # Initial mean/cov
    ukf.x = y0.copy()
    scale = np.maximum(np.abs(y0), 1e-9)

    γ_pH = 5e-5     # ~1% per step on coupled rows
    γ_unc = 0.0
    eps2 = 1e-16

    q_vec = np.zeros(N_STATE)
    for i in range(N_STATE):
        γ = γ_pH if i in coupled_idx else γ_unc
        q_vec[i] = max((γ * scale[i]) ** 2, eps2)
    Q_base = np.diag(q_vec)

    p_vec = np.where(ALGVAR == 1, (0.01 * scale)**2, (0.001 * scale)**2)
    p_vec = np.maximum(p_vec, eps2)
    p_vec[Hfree] = max(p_vec[Hfree], (0.05 * scale[Hfree])**2)
    ukf.P = _make_pos_def(np.diag(p_vec), epsilon=eps2)

    ukf.R = np.atleast_2d(float(r_std) ** 2)

    # Storage
    N = len(t_meas)
    X_hist   = np.zeros((N, N_STATE))
    P_hist   = np.zeros((N, N_STATE, N_STATE))
    pH_pred  = np.zeros(N)
    pH_sigma = np.zeros(N)
    trace_P  = np.zeros(N)
    K_hist   = np.zeros((N, N_STATE))

    for k, (t, z, dt) in enumerate(zip(t_meas, ph_meas, dt_seq)):
        ukf.Q = Q_base * dt
        try:
            ukf.predict(dt=dt)
        except np.linalg.LinAlgError:
            if verbose:
                print(f"[UKF] predict failed at step {k}, repairing P")
            ukf.P = _make_pos_def(ukf.P, epsilon=eps2)
            ukf.predict(dt=dt)

        try:
            ukf.update(z)
        except np.linalg.LinAlgError:
            if verbose:
                print(f"[UKF] update failed at step {k}, repairing P")
            ukf.P = _make_pos_def(ukf.P, epsilon=eps2) + np.eye(N_STATE) * eps2
            ukf.update(z)

        ukf.P = _make_pos_def(ukf.P, epsilon=eps2)

        X_hist[k]   = ukf.x
        P_hist[k]   = ukf.P
        pH_pred[k]  = _hx(ukf.x)[0]
        trace_P[k]  = np.trace(ukf.P)
        K_hist[k]   = ukf.K[:, 0]
        s00 = ukf.S[0, 0]
        pH_sigma[k] = np.sqrt(max(s00, 0.0)) if s00 >= 0 or abs(s00) < 1e-16 else 0.0

    return {"t": t_meas, "x": X_hist, "P": P_hist,
            "pH_pred": pH_pred, "pH_sigma": pH_sigma,
            "trace_P": trace_P, "K": K_hist}

# -------------------------------------------------------------------------
# 10) DATA HELPER
# -------------------------------------------------------------------------
def load_and_sample_ph(
    path: str | Path,
    N_target: int = 100,
    t_final: float = 120.0,
    col: str | None = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """Load pH time-series from CSV/Excel and linearly resample to N_target points over [0, t_final]."""
    path = Path(path)
    if path.suffix.lower() in {".xls", ".xlsx"}:
        df = pd.read_excel(path)
    else:
        df = pd.read_csv(path)
    if col is None:
        col = df.columns[0]
    ph_full = df[col].to_numpy()

    N_full = len(ph_full)
    t_full = np.linspace(0.0, t_final, N_full)
    t_sample = np.linspace(0.0, t_final, N_target)
    pH_sample = np.interp(t_sample, t_full, ph_full)
    return t_sample, pH_sample

# -------------------------------------------------------------------------
# 11)load → UKF → outputs
# -------------------------------------------------------------------------
def run_soft_sensor(
    data_path: str = "data/raw/softsensor_example.csv",
    column: str | None = None,
    t_final: float = 120.0,
    n_samples: int = 121,
    pH0: float = 6.67,
    r_std: float = 0.05,
    run_name: str = "softsensor_run",
    save_pdf: bool = False,
):
    """
    Load pH data, run UKF, save results/figures.

    Returns
    -------
    df : pandas.DataFrame
        Results table (time index, states, pH prediction, uncertainties, etc.)
    outputs : dict
        Paths to CSV and figures directory.
    """
    reports_dir = Path("reports") / run_name
    figures_dir = Path("figures") / run_name
    reports_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # 1) Data
    data_path = Path(data_path)
    if not data_path.exists():
        raise FileNotFoundError(
            f"Data file not found: {data_path}\n"
            "Put your pH series CSV/XLSX at that path or pass --data <path>."
        )
    t, ph_meas = load_and_sample_ph(
        path=data_path, N_target=n_samples, t_final=t_final, col=column
    )

    # 2) Initial state + params
    y0 = make_initials(pH_target=pH0)
    params = default_kin_params()

    # 3) UKF
    est = run_ukf(
        t_meas=t, ph_meas=ph_meas, y0=y0, params=params, r_std=float(r_std)
    )

    # 4) Table (paper-aligned units/labels)
    x_mean = est["x"]
    diag_P = np.diagonal(est["P"], axis1=1, axis2=2)
    x_std  = np.sqrt(diag_P)
    pH_pred  = est["pH_pred"]
    pH_sigma = est["pH_sigma"]
    trace_P  = est["trace_P"]
    K_hist   = est["K"]
    K_norm   = np.linalg.norm(K_hist, axis=1)

    MW_per_nt = 330.0
    MW_RNA = N_TOTAL * MW_per_nt  # g/mol per full transcript

    cols = {
        "time_min": t,
        "pH_meas": ph_meas,
        "pH_pred": pH_pred,
        "pH_pred_minus_2sigma": pH_pred - 2 * pH_sigma,
        "pH_pred_plus_2sigma":  pH_pred + 2 * pH_sigma,
        "trace_P": trace_P,
        "K_norm": K_norm,
    }

    for i, name in enumerate(STATE_NAMES):
        mu = x_mean[:, i].copy()
        sd = x_std[:, i].copy()
        # units for readability
        if name.lower() == "rna":
            mu *= MW_RNA; sd *= MW_RNA        # g/L if original was mol/L
        elif name.endswith("_tot") and "RNA" not in name:
            mu *= 1000.0; sd *= 1000.0        # mol/L -> mM
        cols[name] = mu
        cols[f"{name}_minus_2sigma"] = mu - 2 * sd
        cols[f"{name}_plus_2sigma"]  = mu + 2 * sd
        cols[f"K_{name}"] = K_hist[:, i]

    df = pd.DataFrame(cols).set_index("time_min")
    out_csv = reports_dir / "ivt_ukf_results.csv"
    df.to_csv(out_csv)
    print(f"[OK] results → {out_csv}")

    # 5) Plots
    def _savefig(fig: plt.Figure, base: Path):
        png = base.with_suffix(".png")
        fig.savefig(png, dpi=300, bbox_inches="tight")
        if save_pdf:
            fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)
        print(f"[OK] figure → {png}")

    # pH fit
    fig = plt.figure(figsize=(6, 4))
    ax = fig.gca()
    ax.fill_between(t, df["pH_pred_minus_2sigma"], df["pH_pred_plus_2sigma"], alpha=0.10, label="UKF ±2σ")
    ax.plot(t, df["pH_pred"], lw=2, label="UKF prediction")
    ax.scatter(t, df["pH_meas"], s=12, label="Measurements")
    ax.set_xlabel("Time / min"); ax.set_ylabel("pH"); ax.legend(loc="best")
    _savefig(fig, figures_dir / "pH_fit")

    # diagnostics
    for series, label in [(df["trace_P"], "trace_P"), (df["K_norm"], "K_norm")]:
        fig = plt.figure(figsize=(6, 4))
        ax = fig.gca(); ax.plot(t, series, lw=2)
        ax.set_xlabel("Time / min"); ax.set_ylabel(label)
        _savefig(fig, figures_dir / f"diagnostics_{label}")

    # representative hidden states
    for name in ["RNA", "ATP_tot", "CTP_tot", "GTP_tot", "UTP_tot", "Mg_free"]:
        fig = plt.figure(figsize=(6, 4))
        ax = fig.gca()
        ax.fill_between(t, df[f"{name}_minus_2sigma"], df[f"{name}_plus_2sigma"], alpha=0.20, label=f"{name} ±2σ")
        ax.plot(t, df[name], lw=2, label=name)
        ax.set_xlabel("Time / min"); ax.set_ylabel(name); ax.legend(loc="best")
        _savefig(fig, figures_dir / f"state_{name}")

    print("[DONE] IVT UKF soft sensor run complete.")
    return df, {"csv": out_csv, "fig_dir": figures_dir}

# -------------------------------------------------------------------------
# 12) MAIN (run with defaults)
# -------------------------------------------------------------------------
if __name__ == "__main__":
    # Edit these defaults or call run_soft_sensor(...) from a notebook/script.
    run_soft_sensor(
        data_path="data/softsensor.csv",
        run_name="softsensor_run",
        save_pdf=False,
    )
