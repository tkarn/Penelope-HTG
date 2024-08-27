# Penelope-HTG


## Supplementary Resource to:  

Denkert C et al. 2024 

Dynamics of molecular heterogeneity in high-risk luminal breast cancer - from intrinsic to adaptive subtyping.

************************************************************

## This resource contains the following data regarding the analyses described in the paper:

### Gene signature validations and comparisons based on the HTG gene panel

1. [*HTG-validation_RS-GGI-SET-Pen.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/HTG-validation_RS-GGI-SET-Pen.pdf):  An R Markdown file containing the analyses and figures comparing Recurrence-Score, GGI, SET-Index in TCGA samples based on all genes and the subset of genes available among the 2549 genes from the HTG EdgeSeq Oncology Biomarker Panel measured on the HTG-Molecular EdgeSeq platform.
2. [*HTG-validation_RS-GGI-SET-Pen.R*](https://github.com/tkarn/Penelope-HTG/blob/main/HTG-validation_RS-GGI-SET-Pen.R):  The R-script that generates this R-Markdown file.

#### Input files:
[*HTG-OncBiomarkerPanel_n2559-Genelist.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/HTG-OncBiomarkerPanel_n2559-Genelist.txt): list of genes from HTG-panel.

[*SET-ERPR-genes_PMID_31231679.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/SET-ERPR-genes_PMID_31231679.txt): SET-ER/PR index genelist from PMID 31231679 (Sinn 2019) Suppl.Table 2.

[*Penelope_n355genes_info.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Penelope_n355genes_info.txt): 355 Penelope signature infos.

#### Output files:
[*1_Oncotype-HTG-scatter.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/1_Oncotype-HTG-scatter.svg)
[*2a_GGI-HTG-scatter.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/2a_GGI-HTG-scatter.svg)
[*2b_GGI-HTG-scatter_ERpos.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/2b_GGI-HTG-scatter_ERpos.svg)
[*3a_SET-HTG-scatter.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/3a_SET-HTG-scatter.svg)
[*4_Pen-scatters.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/4_Pen-scatters.svg)



### AIMS subtyping based on the HTG gene panel:

3. [*AIMS-HTG-validation.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/AIMS-HTG-validation.pdf):  An R Markdown file calculating the AIMS subtypes for TCGA samples based on either all genes from the AIMS rules or only on the subset of genes available among the 2549 genes from the HTG EdgeSeq Oncology Biomarker Panel measured on the HTG-Molecular EdgeSeq platform.
4. [*AIMS-HTG-validation.R*](https://github.com/tkarn/Penelope-HTG/blob/main/AIMS-HTG-validation.R):  The R-script that generates this R-Markdown file.




### Relationship of gene clusters and survival:
5. [*2024-08-26-HRscat_BiopResec_335genes_Diff-Colors.pdf*](https://github.com/tkarn/Penelope-HTG/blob/main/2024-08-26-HRscat_BiopResec_335genes_Diff-Colors.pdf):  An R Markdown file of the analysis comparing the prognostic value of individual genes and their association with gene clusters from unsupervised clustering of Penelope-B samples.
6. [*2024-08-26-HRscat_BiopResec_335genes_Diff-Colors.R*](https://github.com/tkarn/Penelope-HTG/blob/main/2024-08-26-HRscat_BiopResec_335genes_Diff-Colors.R):  The R-script that generates this R-Markdown file.

#### Input files:
[*Penelope_n355genes_info.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/Penelope_n355genes_info.txt): 355 Penelope signature infos.

[*HTG-Pathways.xlsx*](https://github.com/tkarn/Penelope-HTG/blob/main/HTG-Pathways.xlsx): Pathway information (HTG and Hallmark) for all genes in the HTG EdgeSeq Oncology Biomarker Panel.

[*HR_quadrant_plot_biopsy_resec_ALL335genes.txt*](https://github.com/tkarn/Penelope-HTG/blob/main/HR_quadrant_plot_biopsy_resec_ALL335genes.txt): Hazard ratios for DFS for all 335 genes from the Penelope-B clusters, based on pre-treatment and post-treatment samples.

#### Output files:
[*1_HRscat_BiopResec_335genes_col-mainclust_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/1_HRscat_BiopResec_335genes_col-mainclust_labs.svg)

[*2_HRscat_BiopResec_335genes_col-mainclust_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/2_HRscat_BiopResec_335genes_col-mainclust_Nolabs.svg)

[*3_HRscat_BiopResec_335genes_col-GeneClass_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/3_HRscat_BiopResec_335genes_col-GeneClass_labs.svg)

[*4_HRscat_BiopResec_335genes_col-GeneClass_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/4_HRscat_BiopResec_335genes_col-GeneClass_Nolabs.svg)

[*5_HRscat_BiopResec_335genes_col-subclust_labs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/5_HRscat_BiopResec_335genes_col-subclust_labs.svg)

[*6_HRscat_BiopResec_335genes_col-subclust_Nolabs.svg*](https://github.com/tkarn/Penelope-HTG/blob/main/6_HRscat_BiopResec_335genes_col-subclust_Nolabs.svg)



************************************************************
