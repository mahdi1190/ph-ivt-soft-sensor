# code/

- `hh_model.py` (Original version found in `notebooks/HH_model_legacy.ipynb` produced by **Shady Hamed**. This version is modified by **Mahdi Ahmed** to be more consistent with the rest of the code.)
- `soft_sensor_main.py` (Produced by **Mahdi Ahmed**. This is the main soft-sensor code used in the paper and should be the primary reference for reproducing the results.)
- `ivt_model.py` (Produced by **Mahdi Ahmed**. 36-state IVT DAE model + UKF routines. Exposes `default_kin_params`, `make_initials`, `run_ivt`, `run_ukf`, and data helpers.)
- `plot_rna_quad.py` (Produced by **Mahdi Ahmed**. Generates a 2×2 RNA yield (g/L) figure: experiment vs UKF vs H–H. Reads `data/all_data_processed.xlsx` and `reports/softsensor_run/*.csv`; optional `data/HH.xlsx`/`data/H_H_NTP.xlsx`.)
- `plot_ph_quad.py` (Produced by **Mahdi Ahmed**. Generates a 2×2 pH figure: experiment vs UKF. Uses per-sheet pH columns from `data/all_data_processed.xlsx` and model CSVs in `reports/softsensor_run/`.)
- `plot_ntp_quad.py` (Produced by **Mahdi Ahmed**. Generates ATP/GTP/CTP/UTP panels (quad or 8-panel): experiment vs UKF vs H–H. Uses `data/all_data_processed.xlsx`, optional H–H workbooks, and `reports/softsensor_run/*.csv`.)
- `plot_species_grid.py` (Produced by **Mahdi Ahmed**. Multi-panel species plots (e.g., PPi, Pi, MgATP/CTP/GTP/UTP) from UKF outputs; saves to `figures/`.)
- `plot_all_pH.py` (Produced by **Mahdi Ahmed**. Legacy/alternate pH plotting entry point; superseded by `plot_ph_quad.py` but kept for completeness.)
- `__init__.py` (Marks this directory as a Python package to allow imports like `from code.ivt_model import run_ukf`.)
