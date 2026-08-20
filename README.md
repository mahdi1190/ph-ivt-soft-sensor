# Soft sensor based on pH for real-time monitoring of mRNA medicines production — Data & Code

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19629299.svg)](https://doi.org/10.5281/zenodo.19629299)

This repository contains the **data**, **analysis code**, and **reproduction instructions** for the manuscript:

> **A soft sensor based on pH for real-time monitoring of mRNA medicine production**  
> Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis. *Digital Discovery*, 2026, **5**, 2899–2914. DOI: [10.1039/D5DD00417A](https://doi.org/10.1039/D5DD00417A)

Archived code & data snapshot: [10.5281/zenodo.19629299](https://doi.org/10.5281/zenodo.19629299) (concept DOI — always resolves to the latest version).

**Affiliations**
- School of Chemical, Materials and Biological Engineering, University of Sheffield, Sheffield S1 3JD, UK.
- Department of Chemical Engineering, Imperial College London, London SW7 2AZ, UK.
  
*First Author:* Mahdi Ahmed (smhamed3@sheffield.ac.uk)
**Corresponding author:** Zoltan Kis (z.kis@sheffield.ac.uk)

## Quick start

1. Clone the repository and create the environment:
   ```bash
   git clone https://github.com/mahdi1190/ph-ivt-soft-sensor.git
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

> Mahdi Ahmed *et al.* “A soft sensor based on pH for real-time monitoring of mRNA medicine production.” *Digital Discovery*, 2026, **5**, 2899–2914. DOI: [10.1039/D5DD00417A](https://doi.org/10.1039/D5DD00417A)

The **code** is released under the **Academic and Research Use License (ARUL)** (see `LICENSE`).  
It is free to use, modify, and redistribute **for non-commercial research and educational purposes only**.  
**Commercial use is not permitted** without a separate, written license agreement.

For commercial licensing enquiries, please contact:
- Commercialisation Team at the University of Sheffield: https://sheffield.ac.uk/commercialisation
- Legal team at the University of Sheffield: ri-contracts@sheffield.ac.uk
- Dr Zoltan Kis, University of Sheffield: z.kis@sheffield.ac.uk

If you use this code in academic work, please **cite the publication** listed above.

### Data licensing

Unless otherwise noted, **data** in `data/` are licensed **CC BY 4.0** (see `DATA_LICENSE.md`).  
Third-party data retain their original licenses and are documented in `DATA.md`.


## How to cite

When citing the **article** and **this repository**, please include both:

- Article:  
  Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis, **A soft sensor based on pH for real-time monitoring of mRNA medicine production**, *Digital Discovery*, 2026, **5**, 2899–2914. DOI: [10.1039/D5DD00417A](https://doi.org/10.1039/D5DD00417A).

- Code & data:  
  Mahdi Ahmed *et al.*, **Soft sensor based on pH for real-time monitoring of mRNA medicines production (data and code)**, Zenodo, 2026, DOI: [10.5281/zenodo.19629299](https://doi.org/10.5281/zenodo.19629299). GitHub: https://github.com/mahdi1190/ph-ivt-soft-sensor.

BibTeX and other formats are available via the repository's [`CITATION.cff`](CITATION.cff) and the [Zenodo record](https://doi.org/10.5281/zenodo.19629299).

## Contact

For questions or access to any restricted data (if applicable), please open an issue or email the corresponding author at z.kis@sheffield.ac.uk.
