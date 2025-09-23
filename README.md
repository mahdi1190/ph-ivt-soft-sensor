# Soft sensor based on pH for real-time monitoring of mRNA medicines production — Data & Code

This repository contains the **data**, **analysis code**, and **reproduction instructions** for the manuscript:

> **Soft sensor based on pH for real-time monitoring of mRNA medicines production**  
> Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis. (*DOI to be added when available.*)

**Affiliations**
- School of Chemical, Materials and Biological Engineering, University of Sheffield, Sheffield S1 3JD, UK.
- Department of Chemical Engineering, Imperial College London, London SW7 2AZ, UK.

**Corresponding author:** Zoltan Kis (z.kis@sheffield.ac.uk)

## Quick start

1. Clone the repository and create the environment:
   ```bash
   git clone https://github.com/mahdi1190/ph-ivt-soft-sensor-DA.git
   cd <repo>
   conda env create -f environment.yml
   conda activate dd-repro
   ```
2. Reproduce the main results:
   ```bash
   jupyter lab
   ```

## Repository layout

```
code/                  # Python modules, Pyomo models, utilities
notebooks/             # Jupyter notebooks to reproduce figures/tables
data/
  raw/                 # Unmodified source data (usually from Zenodo)
  processed/           # Data produced by scripts in code/
  external/            # Third-party reference data (document licenses)
figures/               # Generated figures for the paper
reports/               # Tables and other artefacts used in the manuscript
.github/workflows/     # Optional CI for tests or linting
```

## Data inventory and provenance

See [`DATA.md`](DATA.md) for a table listing each dataset (name, description, size, checksum, source/DOI, license, and how to regenerate).

## Licensing

- **Code** is MIT licensed (see `LICENSE`).
- **Data** are licensed CC BY 4.0 unless otherwise stated (see `DATA_LICENSE.md`). Third-party data retain their original licenses, which are recorded in `DATA.md`.

## How to cite

When citing the **article** and **this repository**, please include both:

- Article (*update if accepted*):  
  Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis, **Soft sensor based on pH for real-time monitoring of mRNA medicines production**, 2025, DOI: tba.

- Code & data:  
  Mahdi Ahmed *et al.*, **Soft sensor based on pH for real-time monitoring of mRNA medicines production (data and code)**, GitHub, 2025, https://github.com/mahdi1190/ph-ivt-soft-sensor-DA (DOI: tba).

BibTeX and other formats will be available via the repository's `CITATION.cff` and the Zenodo record.

## Contact

For questions or access to any restricted data (if applicable), please open an issue or email the corresponding author at z.kis@sheffield.ac.uk.
