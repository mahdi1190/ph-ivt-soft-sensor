# Soft sensor based on pH for real-time monitoring of mRNA medicines production — Data & Code

This repository contains the **data**, **analysis code**, and **reproduction instructions** for the manuscript:

> **Soft sensor based on pH for real-time monitoring of mRNA medicines production**  
> Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis. (*DOI to be added when available.*)

**Affiliations**
- School of Chemical, Materials and Biological Engineering, University of Sheffield, Sheffield S1 3JD, UK.
- Department of Chemical Engineering, Imperial College London, London SW7 2AZ, UK.
  
*First Author:* Mahdi Ahmed (smhamed3@sheffield.ac.uk)
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
code/                  # Python code for the models themselves as well as any relevant plotting code
notebooks/             # Jupyter notebooks to reproduce figures/tables and legacy code
data/                  # Unless specified, this will constain processed data used by the code
  raw/                 # Unmodified source data
  external/            # Third-party reference data (document licenses)
figures/               # Generated figures for the paper
reports/               # Tables and other artefacts used in the manuscript
.github/workflows/     # Optional CI for tests or linting
```

## Data inventory and provenance

See [`DATA.md`](DATA.md) for a table listing each dataset (name, description, source/DOI, license).

## License and Usage
![License: ARUL](https://img.shields.io/badge/License-ARUL-blue.svg)
This repository provides reference code supporting the publication:

> Mahdi Ahmed *et al.* “Soft sensor based on pH for real-time monitoring of mRNA medicines production.” [Journal / Preprint], 2025.  
> *(Update with full citation when available.)*

The **code** is released under the **Academic and Research Use License (ARUL)** (see `LICENSE`).  
It is free to use, modify, and redistribute **for non-commercial research and educational purposes only**.  
**Commercial use is not permitted** without a separate, written license agreement.

For commercial licensing enquiries, please contact:
- Commercialisation team at the University of Sheffield
- Legal team at the University of Sheffield: ri-contracts@sheffield.ac.uk
- Dr Zoltan Kis, University of Sheffield: z.kis@sheffield.ac.uk

If you use this code in academic work, please **cite the publication** listed above.

### Data licensing

Unless otherwise noted, **data** in `data/` are licensed **CC BY 4.0** (see `DATA_LICENSE.md`).  
Third-party data retain their original licenses and are documented in `DATA.md`.


## How to cite

When citing the **article** and **this repository**, please include both:

- Article (*update if accepted*):  
  Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis, **Soft sensor based on pH for real-time monitoring of mRNA medicines production**, 2025, DOI: tba.

- Code & data:  
  Mahdi Ahmed *et al.*, **Soft sensor based on pH for real-time monitoring of mRNA medicines production (data and code)**, GitHub, 2025, https://github.com/mahdi1190/ph-ivt-soft-sensor (DOI: tba).

BibTeX and other formats will be available via the repository's `CITATION.cff` and the Zenodo record.

## Contact

For questions or access to any restricted data (if applicable), please open an issue or email the corresponding author at z.kis@sheffield.ac.uk.
