# Proteome-wide quantification of protein turnover in frog and fly embryos reveals divergent strategies of maternal inheritance

Analysis code and data for the paper of the same title (Cruz et al.). Each analysis is provided as a rendered document (`analysis/*.md`, code with its output) and as the source R Markdown (`scripts/*.Rmd`).

## Data

- Mass spectrometry data: ProteomeXchange/PRIDE **PXD081127** (DDA) and **PXD084487** (DIA).
- `raw/Data/` — search outputs (compressed as `.xz`) and derived intermediates.
- `raw/Files/` — fits, reference tables, absolute concentrations, FASTA databases, and the supplementary tables (`Files/Supp_Tables/Supplementary_Tables_S1-S7.xlsx`).
- `raw/Systems/` — InterProScan and IUPred2A annotations.

The repository is about 1.2 GB on disk (~800 MB to clone). To read the analyses without cloning, open the documents in `analysis/`.

## Analyses

| # | Document | Content |
|---|----------|---------|
| 1 | Frog_GB_NYS-Proteomics | *Xenopus* normalization set and yolk-normalized dataset |
| 2 | Frog_AA_Incorporation_Early-T1 | *Xenopus* 2-cell series: amino acid labeling, peptide filtering |
| 3 | Frog_AA_Incorporation_MBT-T2 | *Xenopus* gastrulation series: amino acid labeling, peptide filtering |
| 4 | Fly_AA_Incorporation_MBT-T2 | *Drosophila* gastrulation series: amino acid labeling, peptide filtering |
| 5 | Frog_2C_SimpModel | *Xenopus* 2-cell: turnover model fitting |
| 6 | Frog_Gast_SimpModel | *Xenopus* gastrulation: turnover model fitting |
| 7 | Fly_Gast_SimpModel | *Drosophila* gastrulation: turnover model fitting |
| 8 | turnover_abs_quant | Absolute protein concentrations |
| 9 | Comp_Frog_Early-MBT | *Xenopus* 2-cell vs. gastrulation turnover |
| 10 | AA-Tracking_Fly_and_Frog | Amino acid supply and demand |
| 11 | Turnover_Species_Comparison | Frog vs. fly turnover at single-copy orthologs |
| 12 | fly_turnover_analysis | *Drosophila* turnover: functional and structural annotation |
| 13 | Frog_Early_Analysis | *Xenopus* turnover: functional and structural annotation |
| 14 | turnover_and_absolute | Turnover vs. absolute abundance |

Rendered versions: `analysis/<Document>.md`. Every intermediate file is included, so each document can be run on its own; the order above follows the pipeline.

## Running the scripts

1. Set the working directory to `raw/` (paths in the scripts are relative to it): `knitr::opts_knit$set(root.dir = "path/to/raw")`, or `setwd()` when running interactively.
2. Decompress the raw search outputs first, from inside `raw/`: `find . -name "*.xz" -exec unxz {} \;` (needed only for the amino acid incorporation documents and `Frog_GB_NYS-Proteomics`).

## Notes

- `Frog_GB_NYS-Proteomics` and `Frog_AA_Incorporation_Early-T1` read each other's outputs; the shipped intermediates cover this. To regenerate both from scratch, run the amino acid matrix step of `Frog_AA_Incorporation_Early-T1` first, then `Frog_GB_NYS-Proteomics`, then the rest.
- "MBT" in file and script names refers to the gastrulation time series.
- Search databases are provided under both `raw/Files/Reference/` and `raw/Files/FASTA/`; the contaminant-containing versions in `FASTA/DIA-NN/` were used for the DIA-NN searches.
- The *Xenopus* normalization reference set (Table S1) is `raw/Data/XLA_Norm/XLA-O18_YolkNormSet_T8-NYS-Decay.csv`; `Frog_GB_NYS-Proteomics` writes its candidate list to `..._candidates.csv`.
- *k*-means, the permutation null and fgsea are seeded (`set.seed(123)`, `clusterSetRNGStream(cl, 123)`).

## Software

R 4.5.3 (package versions in `sessionInfo()` at the end of each rendered document); GFY (Harvard University license) for peptide identification and reporter-ion quantification; DIA-NN 2.5.0 for DIA searches; OrthoFinder for the single-copy ortholog table `raw/Data/XLA-Dmel_orthologue_pairs.csv`.
