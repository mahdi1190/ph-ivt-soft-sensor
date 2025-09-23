#!/usr/bin/env python3
"""
all_data_loader.py
===========

Helper script to deal with the nested data structure of all_data_processed.xlsx

Author: Mahdi Ahmed
License: MIT
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple
import numpy as np
import pandas as pd

# Which pH column to use for each sheet/run
PH_PICK: Dict[str, str] = {
    "egfphepes": "ph2",
    "csptris":   "ph3",
    "csphepes":  "ph2",
    "egfptris":  "ph1",
}

REQUIRED = ["time", "rna2", "ATP_tot", "GTP_tot", "CTP_tot", "UTP_tot"]
OPTIONAL = ["temp", "ph1", "ph2", "ph3", "rna1", "rna3", "qubit1", "qubit2"]

def _clean_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df.columns = [c.strip().replace(" ", "_") for c in df.columns] 
    return df

def _load_sheet_like(path: str | Path, sheet: str) -> pd.DataFrame:
    path = Path(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, sheet_name=sheet)
    else:
        df = pd.read_csv(path)
    return _clean_columns(df)

def load_run(all_data_path: str | Path, run_id: str) -> pd.DataFrame:
    """
    Return a standardized DataFrame for one run/sheet:
    columns: time_min, pH, RNA, ATP_tot, GTP_tot, CTP_tot, UTP_tot, temp
    """
    run = run_id.strip().lower()
    if run not in PH_PICK:
        raise ValueError(f"Unknown run_id '{run_id}'. Expected one of {list(PH_PICK)}")

    df = _load_sheet_like(all_data_path, sheet=run)
    missing = [c for c in REQUIRED if c not in df.columns]
    if missing:
        raise ValueError(f"Sheet '{run}' is missing required columns {missing}. Got: {list(df.columns)}")

    ph_col = PH_PICK[run]
    if ph_col not in df.columns:
        # tolerate missing ph3 etc.
        raise ValueError(f"Sheet '{run}' does not contain column '{ph_col}'. Available: {list(df.columns)}")

    out = pd.DataFrame({
        "time_min": df["time"].astype(float),
        "pH":       df[ph_col].astype(float),
        "RNA":      df["rna2"].astype(float),    #rna 2 is used in this case only 
        "ATP_tot":  df["ATP_tot"].astype(float),
        "GTP_tot":  df["GTP_tot"].astype(float),
        "CTP_tot":  df["CTP_tot"].astype(float),
        "UTP_tot":  df["UTP_tot"].astype(float),
    })
    if "temp" in df.columns:
        out["temp"] = df["temp"].astype(float)
    out = out.sort_values("time_min").reset_index(drop=True)
    return out

def build_master_from_all_data(all_data_path: str | Path) -> pd.DataFrame:
    """
    Build a tidy/long master table (experimental only) from the 4 sheets.
    Columns: run_id, source, time_min, variable, value
    """
    rows = []
    for run in ["egfphepes", "egfptris", "csphepes", "csptris"]:
        d = load_run(all_data_path, run)
        for var in ["ATP_tot","GTP_tot","CTP_tot","UTP_tot","pH","RNA"]:
            series = d["pH"] if var == "pH" else d[var if var in d.columns else var]
            for t, v in zip(d["time_min"].to_numpy(), series.to_numpy()):
                rows.append({"run_id": run, "source": "exp", "time_min": t,
                             "variable": ("pH" if var=="pH" else var.replace("_tot","").upper()),
                             "value": float(v)})
    master = pd.DataFrame(rows).sort_values(["run_id","source","time_min","variable"])
    return master
