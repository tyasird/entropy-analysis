# Differential Interactome & Network Entropy Analysis

This repository provides tools to fetch ovarian cancer gene expression datasets, perform quality control and normalization checks, and apply entropy‑based network/interactome analyses to prioritize candidate biomarkers.

## Table of Contents

- [Motivation](#motivation)  
- [Repository Structure](#repository-structure)  
- [Installation](#installation)  
- [Data Fetching](#data-fetching)  
- [Normalization / QC Checks](#normalization--qc-checks)  
- [Entropy Analysis Pipeline](#entropy-analysis-pipeline)  
- [Inputs & Outputs](#inputs--outputs)  
- [Usage Examples](#usage-examples)  
- [Reproducibility Tips](#reproducibility-tips)  
- [Dependencies](#dependencies)  
- [License & Citation](#license--citation)  
- [Contact](#contact)

---

## Motivation

Changes in network or expression heterogeneity between conditions (e.g. tumor vs normal) can reveal dysregulated modules or key driver genes. Using entropy (information‑theoretic) metrics on interactomes or coexpression graphs helps detect these shifts systematically.

---

## Repository Structure

| File / Folder         | Purpose |
|------------------------|---------|
| `fetch_NCBI.R`          | Script to download GEO / NCBI datasets (expression & metadata) |
| `fetch_TCGA.R`          | Script to fetch and process TCGA ovarian cancer (OV) data |
| `functions.py`          | Core Python utilities: I/O, graph construction, entropy calculations, plotting |
| `norm_check.py`         | Perform QC and normalization diagnostics on expression data |
| `yeni_alg.py`           | Main driver for entropy / differential interactome analysis |
| `data/`                  | (Generated) folder to store downloaded and processed datasets |
| `results/`               | (Generated) folder for outputs: rankings, figures, diagnostics |

---

## Installation

### Python

```bash
python3 -m venv .venv
source .venv/bin/activate   # On Windows: .venv\Scripts\activate
pip install --upgrade pip
pip install numpy pandas scipy networkx matplotlib seaborn tqdm
```

If you have a `requirements.txt`, you can also do:

```bash
pip install -r requirements.txt
```

### R

Run within R (or Rscript) to install needed packages:

```r
install.packages(c("tidyverse", "data.table"))
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "GEOquery"))
```

---

## Data Fetching

### From GEO / NCBI

```bash
Rscript fetch_NCBI.R --accession GSE12345 --outdir data/NCBI
```

This should produce expression and metadata files under `data/NCBI/`.

### From TCGA

```bash
Rscript fetch_TCGA.R --project TCGA-OV --datatype HTSeq --outdir data/TCGA
```

This will download expression and clinical data for ovarian cancer, process them, and save them under `data/TCGA/`.

---

## Normalization / QC Checks

Before running entropy analyses, check that your expression data is suitable (e.g. proper distributions, no extreme outliers, consistent library sizes):

```bash
python norm_check.py   --expr data/TCGA/ov_expression.tsv   --meta data/TCGA/clinical.tsv   --out results/qc
```

Expected outputs: density plots, boxplots, library-size histograms, suggestions for log / transformation, sample QC.

---

## Entropy Analysis Pipeline

Run the main algorithm to build differential interactomes and compute entropy changes:

```bash
python yeni_alg.py   --expr_tumor data/TCGA/ov_tumor.tsv   --expr_normal data/TCGA/ov_normal.tsv   --ppi data/resources/ppi.tsv   --out results/entropy_analysis
```

### What the script does:

1. Load expression matrices (tumor & normal) and PPI / interaction network.  
2. Harmonize gene IDs and filter for common genes.  
3. Build a differential interactome: weight edges by coexpression, change between conditions.  
4. Compute entropy metrics (e.g. per gene, module) in each condition.  
5. Compute **Δentropy** or derived scores.  
6. Rank genes/modules by significance or effect size.  
7. Generate output files (rank tables, figures, network visualizations).

---

## Inputs & Outputs

### Inputs

- **Expression matrices**: genes (rows) × samples (columns), in TSV/CSV format.  
- **PPIs / network file**: edge list format, e.g.:

  ```
  geneA    geneB    weight
  geneC    geneD    1.0
  ...
  ```

  Weight optional.  
- **Metadata (optional)**: sample traits / labels (group, clinical variables).

### Outputs

- `results/entropy_analysis/rankings.csv` — list of genes or modules with Δentropy, scores, ranks  
- `results/entropy_analysis/figures/` — plots (entropy distributions, heatmaps, network diagrams)  
- `results/qc/` — QC and normalization diagnostics  

---

## Usage Examples

```bash
# Fetch data
Rscript fetch_TCGA.R --project TCGA-OV --datatype HTSeq --outdir data/TCGA

# QC
python norm_check.py   --expr data/TCGA/ov_expression.tsv   --meta data/TCGA/clinical.tsv   --out results/qc

# Entropy analysis
python yeni_alg.py   --expr_tumor data/TCGA/ov_tumor.tsv   --expr_normal data/TCGA/ov_normal.tsv   --ppi data/resources/ppi.tsv   --out results/entropy_analysis
```

---

## Reproducibility Tips

- Use deterministic seeds for random steps: `numpy.random.seed(...)`  
- Ensure consistent gene identifiers (symbols vs Ensembl IDs)  
- Filter low-expression genes before graph construction  
- Record versions (`pip freeze > requirements.txt`, `sessionInfo()` in R)

---

## Dependencies

- **Python**: `numpy`, `pandas`, `scipy`, `networkx`, `matplotlib`, `seaborn`, `tqdm`  
- **R**: `tidyverse`, `data.table`, `TCGAbiolinks`, `GEOquery`, `SummarizedExperiment`

---

## License & Citation

Add a `LICENSE` file if publishing (e.g. MIT). Cite this repo in publications if used.

---

## Contact

Author: Talip Yasir Demirtas ([@tyasird](https://github.com/tyasird))  

Please open GitHub Issues for bugs or feature requests.
