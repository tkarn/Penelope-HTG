## Penelope-HTG / Fig2

## Supplementary Resource to:  


Denkert C et al. 2025, Cancer Cell 43, 1–16, February 10, 2025 
https://doi.org/10.1016/j.ccell.2025.01.002

Dynamics of molecular heterogeneity in high-risk luminal breast cancer - From intrinsic to adaptive subtyping.

PMID: 39933898

************************************************************

## This subdirectory contains the code and data to generate Figure 2 from the paper:

### (E) Sankey diagram. Paired analysis of 540 tumors before and after therapy showing the transition from LumB to LumA as the main alteration during neoadjuvant therapy. 
### (H) Sankey diagram of longitudinal AIMS-subtypes in n=29 tumors with no response to chemotherapy and three longitudinal tumor samples indicating a transition from LumB to LumA in the chemotherapy phase that is reversed in metastatic disease. 

1. [*sankey_BR_AIMS_Fig2E.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/sankey_BR_AIMS_Fig2E.pdf):  An R Markdown file of the paired analysis of 540 tumors before and after therapy showing the transition from LumB to LumA as the main alteration during neoadjuvant therapy in a Sankey diagram.
2. [*sankey_BR_AIMS_Fig2E.R*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/sankey_BR_AIMS_Fig2E.R):  The R-script that generates this R-Markdown file.
3. [*sankey_BRM_AIMS_Fig2H.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/sankey_BRM_AIMS_Fig2H.pdf):  An R Markdown file of the paired analysis of longitudinal AIMS-subtypes in n=29 tumors with metastatic samples in a Sankey diagram.

#### Input files:
[*BR_AIMSnew.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/BR_AIMSnew.txt): pre-Tx and post-Tx AIMS classification of 540 paired samples.
[*BRM_AIMSnew.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/BRM_AIMSnew.txt): pre-Tx, post-Tx, and metastasis AIMS classification of 29 paired samples.

#### Output files:
[*Fig2E_Sankey_plotBR_AIMS.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/Fig2E_Sankey_plotBR_AIMS.pdf)

[*Fig2H_Sankey_plotBRM_AIMS_rev.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/Fig2/Fig2H_Sankey_plotBRM_AIMS_rev.pdf)

************************************************************
