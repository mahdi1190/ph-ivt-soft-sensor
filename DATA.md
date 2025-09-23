# DATA.md — Dataset inventory & access

**Manuscript:** Soft sensor based on pH for real-time monitoring of mRNA medicines production  
**Authors:** Mahdi Ahmed, Shady Hamed, Ricardo Cardoso, Charley Kenyon, Manoj Pohare, Mabrouka Maamra, Mark Dickman, Joan Cordiner, Zoltan Kis  
**Corresponding author:** Zoltan Kis (z.kis@sheffield.ac.uk)


This table lists *all* datasets used or generated in the paper. For each item, provide a DOI/URL and access conditions.

| ID | Path | Description | Source / DOI | License | Size (MB) | SHA256 | Access | Regenerate command |
|----|------|-------------|--------------|---------|-----------|--------|--------|--------------------|
| D1 | data/raw/<file> | <short desc> | <DOI or URL> | <CC BY 4.0 / proprietary / etc.> | <size> | <checksum> | public | `python code/pipeline.py --stage fetch --item D1` |
| D2 | data/processed/<file> | <derived desc> | generated | CC BY 4.0 | <size> | <checksum> | public | `python code/pipeline.py --stage process --from D1` |