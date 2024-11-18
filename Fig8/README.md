## Penelope-HTG / Fig8

## Supplementary Resource to:  

Denkert C et al. 2024 

Dynamics of molecular heterogeneity in high-risk luminal breast cancer - from intrinsic to adaptive subtyping.

************************************************************

## This subdirectory contains the code and data to generate Figure 8 from the paper:

### Identification of prognostic groups by additional UMAP analysis.

1. [*UMAP_Fig8a.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_Fig8a.pdf):  An R Markdown file generating an UMAP for 1080 tumor samples (540 pre-Tx and 540 paired post-Tx samples) with euclidean metric for the set of 335 differentially regulated genes between pre-Tx and post-TX samples, coloured by AIMS classification.
2. [*UMAP_Fig8a.R*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_Fig8a.R):  The R-script that generates this R-Markdown file.
3. [*UMAP_Fig8b.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_Fig8b.pdf):  An R Markdown file generating the same UMAP coloured by Adaptive Subtypes (AC) classification.
4. [*UMAP_Fig8b.R*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_Fig8b.R):  The R-script that generates this R-Markdown file.

#### Input files:
[*UMAP335_Biop_Res.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP335_Biop_Res.txt): Gene expression data.

[*UMAP_sample_info.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_sample_info.txt): Sample information.

#### Output files:
[*UMAP_AIMS_Fig8A.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_AIMS_Fig8A.svg)

[*UMAP_AC_clusters_Fig8B.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig8/UMAP_AC_clusters_Fig8B.svg)

************************************************************
