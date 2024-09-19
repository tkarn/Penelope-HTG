# HEADER ####
#
# Version: 2024-08-26
#
# Comparison of AIMS results in TCGA-BRCA-RNA-Seq based on 
# complete gene list and the subset available in the HTG-Panel
#
#

#
# SETUP ####
#' 
Sys.setenv(lang = "en_US")

#' *Install required packages if missing* -----------------------------------


# Package names from CRAN
packs <- c("dplyr", "AIMS", "crosstable")

# Install packages not yet installed
installed_packages <- packs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packs[!installed_packages])
}


# Package names from Bioconductor
bcpacks <- c("cBioPortalData")

# Install bc-packages if not yet installed from Bioconductor
installed_packages <- bcpacks %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(bcpacks[!installed_packages])
}


#' *Load required packages* -----------------------------------

invisible(lapply(packs, library, character.only = TRUE))
invisible(lapply(bcpacks, library, character.only = TRUE))



#'\newpage
# FUNCTION Definitions ####
# ******************************************************
#
# Function tcgaRseqEntrez  ####
#
#
# The function tcgaRseqEntrez obtaines RNA-Seq data of a provided list 
#  of entrezGeneId's from TCGA using the cBioPortal access tools
#  and delivers only the data for BRCA samples with ER status available
#   (ERpos and ERneg samples).
#
# THIS IS A SLIGHT CHANGED VERSION OF FUNCTION tcgaRseqGenelist from
#  the file HTG-validation_RS-GGI-SET-Pen.R
#
# We apply the cBioPortalData package to access data from the cBIO Portal
# at www.cbioportal.org
# This will allow to download RNA-Seq data from the TCGA-BRCA cohort.

library(cBioPortalData)
# First we setup some parameters for the cBioportal-access
# Define api
cbio <- cBioPortal()

# Function definition:
#  (entrezList is a vector of entrezGeneId's)

tcgaRseqEntrez <- function (entrezList) {
  # Download BRCA RNA-Seq data for this list of entrezGeneId's from cBioPortal
  #  as a "MultiAssayExperiemnt" brca_rnaseq
  brca_rnaseq <- cBioPortalData(
    api = cbio,
    studyId = "brca_tcga",
    genes = entrezList, by = "entrezGeneId",
    molecularProfileIds = "brca_tcga_rna_seq_v2_mrna"
    )
  # Extract the RNA-Seq data from the MultiAssayExperiment
  tcgaRseqGenelist <- assay(brca_rnaseq[["brca_tcga_rna_seq_v2_mrna"]])
  
  # Extract the phenotype data for TCGA-samples (by patientId)
  pheno <- colData(brca_rnaseq)
  
  # Extract the link-information between
  # the patientId ("primary") and  the RNA-seq-colnames ("colname")
  # from the MultiAssayExperiment
  sample_info <- unique(sampleMap(brca_rnaseq)[,2:3])
  
  
  # Now we use dplyr functions from tidyR to join
  #    the "ER_STATUS_BY_IHC" from pheno
  #    with the "colname" from sample_info
  # by linking the cases using the patientId == primary
  
  pdata <- as.data.frame(pheno) %>% 
    dplyr::select(patientId, ER_STATUS_BY_IHC) %>% 
    left_join(as.data.frame(sample_info), by = join_by(patientId == primary))
  
  # We can now use pdata to select only samples with information
  #   on ER status
  
  pdata.erstat <- pdata %>% filter(ER_STATUS_BY_IHC %in% c("Positive", "Negative"))
  
  tcgaRseqEntrez <- tcgaRseqGenelist[, colnames(tcgaRseqGenelist)
                                       %in% pdata.erstat$colname]
}

#'\newpage
# Function tcgaRseqEntrezERpos  ####
#
## The function tcgaRseqEntrezERpos obtaines RNA-Seq data of a provided list 
#  of entrezGeneId's from TCGA using the cBioPortal access tools.
#  and delivers only the data of ERpos BRCA samples.
#
# THIS IS A SLIGHT CHANGED VERSION OF FUNCTION tcgaRseqGenelistERpos from
#  the file HTG-validation_RS-GGI-SET-Pen.R
#
# We apply the cBioPortalData package to access data from the cBIO Portal
# at www.cbioportal.org
# This will allow to download RNA-Seq data from the TCGA-BRCA cohort.

library(cBioPortalData)
# First we setup some parameters for the cBioportal-access
# Define api
cbio <- cBioPortal()

# Function definition:
#  (genelist is a vector of gene symbols)

tcgaRseqEntrezERpos <- function (entrezList) {
  # Download BRCA RNA-Seq data for this genelist from cBioPortal
  #  as a "MultiAssayExperiemnt" brca_rnaseq
  brca_rnaseq <- cBioPortalData(
    api = cbio,
    studyId = "brca_tcga",
    genes = entrezList, by = "entrezGeneId",
    molecularProfileIds = "brca_tcga_rna_seq_v2_mrna"
  )
  # Extract the RNA-Seq data from the MultiAssayExperiment
  tcgaRseqGenelist <- assay(brca_rnaseq[["brca_tcga_rna_seq_v2_mrna"]])
  
  # Extract the phenotype data for TCGA-samples (by patientId)
  pheno <- colData(brca_rnaseq)
  
  # Extract the link-information between
  # the patientId ("primary") and  the RNA-seq-colnames ("colname")
  # from the MultiAssayExperiment
  sample_info <- unique(sampleMap(brca_rnaseq)[,2:3])
  
  
  # Now we use dplyr functions from tidyR to join
  #    the "ER_STATUS_BY_IHC" from pheno
  #    with the "colname" from sample_info
  # by linking the cases using the patientId == primary
  
  pdata <- as.data.frame(pheno) %>% 
    dplyr::select(patientId, ER_STATUS_BY_IHC) %>% 
    left_join(as.data.frame(sample_info), by = join_by(patientId == primary))
  
  # We can now use pdata to select only ER-positive samples
  
  pdata.erpos <- pdata %>% filter(ER_STATUS_BY_IHC == "Positive")
  
  tcgaRseqEntrezERpos <- tcgaRseqGenelist[, colnames(tcgaRseqGenelist)
                                                 %in% pdata.erpos$colname]
}



#'\newpage
# DATA IMPORT ####

library(dplyr)

#' *Import list of genes from HTG-panel*

htgprobes <- pull(read.table("HTG-OncBiomarkerPanel_n2559-Genelist.txt",
                       header=FALSE, sep=","))


# Mapping of HTG gene names to entrez IDs, AIMS genes only
# 
HTG.AIMS.gene.map <- as.data.frame(matrix(c(
  "ANXA3", "306",
  "APH1B", "83464",
  "AR", "367",
  "ASPM", "259266",
  "BCL2", "596",
  "BIRC5", "332",
  "C1orf106", "55765",
  "CA12", "771",
  "CAV1", "857",
  "CCNB2", "9133",
  "CDC20", "991",
  "CDH3", "1001",
  "CDKN1C", "1028",
  "CDKN3", "1033",
  "CENPF", "1063",
  "CEP55", "55165",
  "CIRBP", "1153",
  "CKS2", "1164",
  "CNIH4", "29097",
  "COL17A1", "1308",
  "CRYAB", "1410",
  "CSTB", "1476",
  "CX3CL1", "6376",
  "DNAJC12", "56521",
  "ERBB2", "2064",
  "ESR1", "2099",
  "FBP1", "2203",
  "FGFR4", "2264",
  "FMO5", "2330",
  "FOXA1", "3169",
  "FOXC1", "2296",
  "GAMT", "2593",
  "GATA3", "2625",
  "GFRA1", "2674",
  "GSN", "2934",
  "GSTP1", "2950",
  "HPN", "3249",
  "HSPA14", "51182",
  "ID4", "3400",
  "IGF1", "3479",
  "IGFBP6", "3489",
  "IRS1", "3667",
  "ITM2A", "9452",
  "KIF2C", "11004",
  "KIT", "3815",
  "KRT14", "3861",
  "KRT17", "3872",
  "KRT18", "3875",
  "KRT5", "3852",
  "LAMA3", "3909",
  "LYN", "4067",
  "MAD2L1", "4085",
  "MAP2K4", "6416",
  "MAPT", "4137",
  "MCM2", "4171",
  "MELK", "9833",
  "MLPH", "79083",
  "MMP7", "4316",
  "MNAT1", "4331",
  "NAT1", "9",
  "NDC80", "10403",
  "NEK2", "4751",
  "NQO1", "1728",
  "PARP1", "142",
  "PCNA", "5111",
  "PPAP2B", "8613",
  "PRC1", "9055",
  "PRKX", "5613",
  "PTN", "5764",
  "PTTG1", "9232",
  "RACGAP1", "29127",
  "RBBP8", "5932",
  "RFC4", "5984",
  "RRM2", "6241",
  "S100A8", "6279",
  "SCUBE2", "57758",
  "SERPINA3", "12",
  "SFRP1", "6422",
  "SHC2", "25759",
  "SLC39A6", "25800",
  "SPDEF", "25803",
  "STC2", "8614",
  "TFF3", "7033",
  "TK1", "7083",
  "TNFRSF21", "27242",
  "TOP2A", "7153",
  "TSPAN13", "27075",
  "TSPAN7", "7102",
  "TTK", "7272",
  "TYMS", "7298",
  "UBE2C", "11065"),
  ncol = 2L, byrow = T, dimnames = list(NULL, c("HTGname", "entrezID"))), stringsAsFactors = F)

stopifnot(!duplicated(HTG.AIMS.gene.map$HTGname))
stopifnot(!duplicated(HTG.AIMS.gene.map$entrezID))


library(AIMS)

# internal check: entrez IDs used by AIMS
entrezIDs.AIMS <- unique(unlist(strsplit(AIMSmodel$all.pairs, "<", T)))
stopifnot(length(entrezIDs.AIMS) == 151L)
stopifnot(HTG.AIMS.gene.map$entrezID %in% entrezIDs.AIMS)

entrezIDs.AIMS.HTG <- unique(HTG.AIMS.gene.map$entrezID)
# 91 entrez IDs



#'\newpage
# ANALYSIS ####


# AIMS using all entrezIDs for TCGA samples####

# Load RNA-Seq data for all entrezIDs from AIMS package

entrezList <- entrezIDs.AIMS

# RNAseq for all TCGA-BRCA (including ERneg)
tcga.AIMS.Rseq <- tcgaRseqEntrez(entrezList)

# RNAseq for ERpos TCGA-BRCA
tcga.AIMS.Rseq.ERpos <- tcgaRseqEntrezERpos(entrezList)


# ***************************************************************
#  Calculate AIMS groups including ALL 151 entrezIDs from TCGA  *
# ***************************************************************

# All TCGA samples
tcga.AIMS <- applyAIMS(tcga.AIMS.Rseq, rownames(tcga.AIMS.Rseq))

# ERpos subset
tcga.AIMS.ERpos <- applyAIMS(tcga.AIMS.Rseq.ERpos, rownames(tcga.AIMS.Rseq.ERpos))


# ***************************************************************
#  Calculate AIMS groups only for entrezIDs with HTG-data       *
# ***************************************************************

# All TCGA samples
tcga.AIMS.Rseq.HTG <- tcga.AIMS.Rseq[entrezIDs.AIMS.HTG,]

tcga.AIMS.HTG <- applyAIMS(tcga.AIMS.Rseq.HTG, rownames(tcga.AIMS.Rseq.HTG))


# ERpos subset
tcga.AIMS.ERpos.Rseq.HTG <- tcga.AIMS.Rseq.ERpos[entrezIDs.AIMS.HTG,]

tcga.AIMS.ERpos.HTG <- applyAIMS(tcga.AIMS.ERpos.Rseq.HTG,
                                 rownames(tcga.AIMS.ERpos.Rseq.HTG))



# *************************************************************************

#'\newpage
# COMPARE OBTAINED RESULTS FOR AIMS ####

# Finally we compare the results obtained from the entrezID set rules


# All TCGA samples

# All rules
table(tcga.AIMS$cl)
# HTG gene rules only
table(tcga.AIMS.HTG$cl)

# cross table
table(tcga.AIMS$cl, tcga.AIMS.HTG$cl)



# ERpos samples

# All rules
table(tcga.AIMS.ERpos$cl)
# HTG gene rules only
table(tcga.AIMS.ERpos.HTG$cl)

# cross table
table(tcga.AIMS.ERpos$cl, tcga.AIMS.ERpos.HTG$cl)

# *************************************************************************
#'\newpage
# CROSSTAB with nice formatting and percentages

ct.labels = read.table(header=TRUE, text="
  name        label
  aims.all    'Subtype all genes'
  aims.all.pr 'Prob all genes'
  aims.htg    'Subtype HTG'
  aims.htg.pr 'Prob HTG'
  ")



# For all TCGA samples:
#
# Define dataframe with colnames
tcga.AIMS.ct <- as.data.frame(cbind(tcga.AIMS$cl, tcga.AIMS$prob,
                                    tcga.AIMS.HTG$cl, tcga.AIMS.HTG$prob))
colnames(tcga.AIMS.ct) = c("aims.all", "aims.all.pr", "aims.htg", "aims.htg.pr")
# Import labels in tibble
tcga.AIMS.ct <- tcga.AIMS.ct %>% 
  import_labels(ct.labels, name_from="name", label_from="label") %>% 
  as_tibble()


#'\newpage
#' **Cross table with column-percentages (subtype all genes) for all TCGA samples**
crosstable(tcga.AIMS.ct, aims.all, by=aims.htg, total="both", percent_digits=1,
           percent_pattern="{n}\n{p_col}")  %>% 
  as_flextable(compact = TRUE)

#'\newpage
#' **Cross table with row-percentages (subtype htg) for all TCGA samples**
crosstable(tcga.AIMS.ct, aims.all, by=aims.htg,
           total="both",
           percent_digits=1,
           percent_pattern="{n}\n{p_row}")  %>% 
  as_flextable(compact = TRUE)



#'\newpage
# *************************************************************************
# For ERpos subset of samples:

# Define dataframe with colnames
tcga.AIMS.ERpos.ct <- as.data.frame(cbind(tcga.AIMS.ERpos$cl, tcga.AIMS.ERpos$prob, 
                                 tcga.AIMS.ERpos.HTG$cl, tcga.AIMS.ERpos.HTG$prob))
colnames(tcga.AIMS.ERpos.ct) = c("aims.all","aims.all.pr","aims.htg","aims.htg.pr")
# Import labels in tibble
tcga.AIMS.ERpos.ct <- tcga.AIMS.ERpos.ct %>% 
  import_labels(ct.labels, name_from="name", label_from="label") %>% 
  as_tibble()



#' **Cross table with column-percentages (subtype all genes) for ERpos TCGA subset:**
crosstable(tcga.AIMS.ERpos.ct, aims.all, by=aims.htg, total="both",
           percent_digits=1, percent_pattern="{n}\n{p_col}")  %>% 
  as_flextable(compact = TRUE)


#'\newpage
#' **Cross table with row-percentages (subtype htg) for ERpos TCGA subset:**
crosstable(tcga.AIMS.ERpos.ct, aims.all, by=aims.htg, total="both", 
           percent_digits=1, percent_pattern="{n}\n{p_row}")  %>% 
  as_flextable(compact = TRUE)



#'\newpage
# *************************************************************************
# SESSION INFO ####
sessionInfo()

