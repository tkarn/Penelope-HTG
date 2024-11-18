# HEADER ####
#
# Version: 2024-11-18
#
# Figure 2H: Sankey plots of AIMS subtypes in pre-Tx / post-Tx / metastasis
# 
# !! IMPORTANT !!:
#
# The devtools() package must already be installed 
#   for the installation of ggsankey()
#
#
# SETUP ####
#'  
Sys.setenv(lang = "en_US")


#' *Install required packages if missing* -----------------------------------


# Package ggsankey from GitHub using devtools
if ("ggsankey" %in% rownames(installed.packages) == FALSE) {
  devtools::install_github("davidsjoberg/ggsankey")
} 


# Package names for install from CRAN
packs <- c("ggplot2", "dplyr", "ggalluvial", "networkD3")

# Install packages not yet installed
installed_packages <- packs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packs[!installed_packages])
}

#' *Load required packages* -----------------------------------

invisible(library(ggsankey))

invisible(lapply(packs, library, character.only = TRUE))

set.seed(321)

#'\newpage
# IMPORT ####

# Import AIMS subtypes of pre-Tx biopsy, post-Tx resect, and metast. sample

sank_df <- read.delim("BRM_AIMS.txt", header=T, stringsAsFactors = T, skipNul=T)

# Analysis / Plots ####

df <- sank_df %>%
  make_long(pre.Tx,post.Tx, metastasis) |> mutate(next_node = forcats::fct_inorder(next_node))

# counts and percentages
  TotalCount = nrow(sank_df)
  dagg <- df%>%
    dplyr::group_by(node)%>%
    tally()
  
  dagg <- dagg%>%
    dplyr::group_by(node)%>%
    dplyr::mutate(pct = n/TotalCount) 


# visualising plot
  df2 <- merge(df, dagg,   by.x = 'node',by.y = 'node')
  pl <- ggplot(df2, aes(x = x
                        , node = node
                        , next_x = next_x
                        , next_node = next_node
                        , fill = factor(node),
                        , label = paste0(node," n=", n, ' (',  round(pct* 100,1), '%)' )) )
  
  pl <- pl +geom_sankey(width = 1/4, flow.alpha = 1.5,  node.color = "gray40", show.legend = TRUE)+
        geom_sankey_label(size = 3, color = "black", fill= "white", hjust = -0.2)+
        theme_bw()+
        theme(legend.position = "none")+
        theme(axis.title = element_blank(),
                           axis.text.y = element_blank(),
                           axis.ticks = element_blank(),
                           panel.grid = element_blank())+
        
        labs(fill = 'Nodes')+
        labs(title = "Sankey diagram")+
        scale_fill_manual(values =c('b.Basal'="red",'b.HER2E'="magenta",'b.LumA'='darkblue','b.LumB'="skyblue",'b.NormL'="green",
                                               'r.Basal'="red",'r.HER2E'="magenta",'r.LumA'='darkblue','r.LumB'="skyblue",'r.NormL'="green",
                                               'm.Basal'="red",'m.HER2E'="magenta",'m.LumA'='darkblue','m.LumB'="skyblue",'mr.NormL'="green"))
  pl
ggsave( "./Fig2H_Sankey_plotBRM_AIMS_rev.pdf", pointsize = 12, bg = "white")



# SESSION INFO ####
sessionInfo()