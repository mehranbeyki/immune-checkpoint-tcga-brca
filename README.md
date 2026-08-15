# Single-Cell RNA-seq Analysis of Peripheral Blood in HR+ Breast Cancer Patients Treated with Immunotherapy

> **Status: Work in progress.** This repository is under active development. Notebooks, scripts, and results are being iterated on and may change.

## Overview

This project analyzes a public single-cell RNA sequencing (scRNA-seq) dataset investigating peripheral immune dynamics in hormone receptor-positive (HR+) breast cancer patients treated with neoadjuvant **nab-paclitaxel + pembrolizumab**. The data combine **CITE-seq** (gene expression + hashtag/antibody capture) and **TCR-seq** (V(D)J) libraries generated on the 10x Genomics Chromium platform from longitudinal peripheral blood mononuclear cell (PBMC) samples across multiple patients and treatment timepoints.

The overall aim is to reproduce and extend the published analysis by performing sample demultiplexing, quality control, integration, clustering, and cell type annotation, in order to explore peripheral blood correlates of immunotherapy response.

- **GEO Accession:** [GSE300475](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300475)
- **BioProject:** [PRJNA1281033](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1281033)
- **Organism:** *Homo sapiens*
- **Platform:** Illumina NovaSeq 6000 (GPL24676)
- **Assay types:** 5' scRNA-seq, Feature Barcode (HTO), V(D)J (TCR)
- **Samples:** 32 GEO records covering multiple patients (PT1, PT5, PT6, PT7, PT11, PT13, PT15) across timepoints (week1, week3, week7, week15), each hashtag-multiplexed (CITE-seq HTO)
- **Original study contact:** Balko Lab, Vanderbilt University
- **Publication:** [PubMed 40610460](https://www.ncbi.nlm.nih.gov/pubmed/40610460)

## Study Summary

The limited clinical benefit of immune checkpoint inhibitors in breast cancer highlights the need for predictive biomarkers that can minimize risk and maximize benefit for patients. The original study performed single-cell RNA and TCR sequencing on PBMCs to monitor peripheral immune dynamics in an exploratory cohort of HR+ breast cancer patients undergoing neoadjuvant immunotherapy, aiming to identify candidate peripheral blood biomarkers of treatment response.

## Data & Results Versioning

Raw data, intermediate files, figures, and result tables are **not stored in this Git repository**. They are version-controlled with **DVC** and hosted on DagsHub:

🔗 **https://dagshub.com/mehranbeyki/scrna-immune-checkpoint-response**

To pull the tracked data/results locally, clone the repository and run the standard DVC pull workflow (`dvc pull`) after configuring the DagsHub remote — see the DagsHub project page for setup instructions.

The original raw sequencing/count files are also available directly from GEO:

- Supplementary raw data: `GSE300475_RAW.tar` (per-sample CSV/MTX/TSV files, 10x Genomics format)
- Feature reference: `GSE300475_feature_ref.xlsx`
- Raw FASTQ reads: via [SRA Run Selector](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1281033)

## Project Structure

The structure below is an **illustrative example** based on the current notebooks/scripts — the actual layout will be finalized as the project progresses.

```
.
├── data/
│   ├── raw/
│   │   ├── GSE300475_RAW/        # Original files extracted from GSE300475_RAW.tar
│   │   ├── gex/                  # Per-sample gene expression + HTO folders (organized by file_org.py)
│   │   │   └── scRNAseq PT1/     #   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
│   │   └── vdj/                  # Per-sample V(D)J contig annotation files
│   ├── mapping_table.csv         # GSM accession → sample title mapping (from GEOquery metadata)
│   └── processed/                # Intermediate & final Seurat objects (.rds), tracked via DVC
├── notebooks/
│   ├── 01_data_download_preprocessing.ipynb   # Download, HTO demultiplexing, per-sample Seurat objects
│   ├── 02_qc_clustering_R.ipynb               # QC, patient/timepoint metadata mapping, SCTransform, Harmony integration, clustering
│   └── 03_celltype_annotation_.ipynb          # Marker-based cell type annotation
├── scripts/
│   ├── file_org.py                # Organizes raw GEO files into per-sample 10x-style folders (gex/, vdj/)
│   └── check_gzip_features.py     # Diagnostic script to inspect gzipped features.tsv.gz files
├── results/
│   ├── figures/                   # QC plots, UMAPs, dot plots, etc. (tracked via DVC)
│   └── tables/                    # Marker gene tables, QC summary tables (tracked via DVC)
├── README.md
└── dvc.yaml / .dvc files          # DVC pipeline and data tracking configuration
```

## Workflow

The pipeline combines a Python preprocessing step with an R-based analysis pipeline (run as Jupyter notebooks with an R kernel).

### 1. Raw data organization (Python)

After downloading and extracting `GSE300475_RAW.tar` from GEO, `scripts/file_org.py` reorganizes the flat archive of per-GSM files into per-sample folders following the 10x Genomics convention (`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` for gene expression samples, and `<sample>_all_contig_annotations.csv.gz` for V(D)J samples). Sample folder names are resolved via `data/mapping_table.csv` (GSM accession → sample title), which is generated from GEO metadata in notebook `01`.

```bash
python scripts/file_org.py
```

`scripts/check_gzip_features.py` is a small diagnostic utility used to inspect the contents of a sample's `features.tsv.gz` file (e.g., to confirm feature counts/formatting) before loading it into R.

### 2. Data download, metadata mapping & demultiplexing — `01_data_download_preprocessing.ipynb` (R)

- Downloads GEO supplementary files and series metadata via `GEOquery` (`getGEOSuppFiles`, `getGEO`)
- Builds `mapping_table.csv` linking GEO accessions, sample titles, and source names
- Loads each sample with `Seurat::Read10X`, creates a `Seurat` object (gene expression) with an additional `HTO` assay (antibody capture)
- Normalizes HTO counts with CLR (Centered Log-Ratio) and demultiplexes hashtags with `HTODemux`
- Visualizes HTO signal (ridge plots, heatmaps, violin plots) and retains **singlets** only
- Loops the full pipeline across all sample folders and merges all per-sample singlet objects into a single Seurat object (`merged_singlets.rds`)

### 3. QC, patient/timepoint mapping, integration & clustering — `02_qc_clustering_R.ipynb` (R)

- Computes mitochondrial content (`percent.mt`) and inspects QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`) by sample
- Builds a manual **hashtag → patient/timepoint mapping table** and joins it to cell-level metadata to annotate each cell with `patient_id` and `timepoint`
- Flags the PT15 replicate run (`sequencing_run` column) and performs a QC comparison between the original and replicate runs
- Applies QC filtering (`200 < nFeature_RNA < 6000`, `nCount_RNA < 30000`, `percent.mt < 15`)
- Runs `SCTransform` per patient (regressing out `percent.mt`), then merges and selects integration features
- Integrates samples with **Harmony** (`IntegrateLayers`), computes PCA/UMAP, and clusters cells (`FindNeighbors` + `FindClusters`)
- Saves diagnostic UMAPs (by patient, sequencing run, cluster, timepoint) and the final integrated/clustered object (`merged_integrated_clustered.rds`)

### 4. Cell type annotation — `03_celltype_annotation_.ipynb` (R)

- Loads the integrated/clustered Seurat object
- Scores clusters against a curated panel of PBMC lineage markers (T cell, CD4/CD8, Treg, exhaustion markers, NK, monocyte, B cell, plasma cell, dendritic cell, platelet, red blood cell) via `DotPlot`
- Identifies cluster marker genes with `FindAllMarkers` (SCT assay) and exports full and top-5-per-cluster marker tables
- Saves the annotated/prepared object for downstream analysis

## Requirements

**R** (primary analysis environment, R ≥ 4.5)

Key packages used across the notebooks:
`Seurat`, `SingleR`, `celldex`, `GEOquery`, `dplyr`, `ggplot2`, `patchwork`, `RColorBrewer`, `glmGamPoi`, `harmony`, `future`

**Python** (used for raw data organization only)

`pandas`, plus standard library modules `os`, `re`, `shutil`, `gzip`

**Data/results versioning**

`dvc` (with a DagsHub remote configured — see [Data & Results Versioning](#data--results-versioning))

## Usage

1. Clone this repository.
2. Configure the DVC remote and pull tracked data/results (see the [DagsHub project page](https://dagshub.com/mehranbeyki/scrna-immune-checkpoint-response)), **or** download the raw files directly from GEO and place them under `data/raw/`.
3. Run `scripts/file_org.py` to organize the raw GEO files into per-sample 10x-style folders.
4. Run the notebooks in order:
   1. `notebooks/01_data_download_preprocessing.ipynb`
   2. `notebooks/02_qc_clustering_R.ipynb`
   3. `notebooks/03_celltype_annotation_.ipynb`

## Citation

If you use this dataset, please cite the associated publication:

- PubMed: [40610460](https://www.ncbi.nlm.nih.gov/pubmed/40610460)

Please also cite the GEO series:

> Sun X, Axelrod ML, et al. Single cell RNA sequencing for longitudinal human peripheral blood from HR+ breast cancer patients treated with immunotherapy. GEO accession GSE300475.

## Contact

For questions about the original dataset, contact the submitting lab:

- **Justin Balko** — Vanderbilt University — justin.balko@vumc.org

For questions about this analysis repository, please open an issue on GitHub.
