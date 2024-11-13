# HEADER ####
#
# Version: 2024-11-13
#
# Figure S5: 
#   Analysis of prognostic gene sets in pre-Tx biopsies 
#   and post-Tx residual tumors.
# 
#
#
#
# SETUP ####
#'  
Sys.setenv(lang = "en_US")

#' *Install required packages if missing* -----------------------------------

# Package names from CRAN
packs <- c("tidyverse", "ggrepel")

# Install packages not yet installed
installed_packages <- packs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packs[!installed_packages])
}

#' *Load required packages* -----------------------------------

invisible(lapply(packs, library, character.only = TRUE))



#'\newpage
# FUNCTION Definitions ####
# ******************************************************
#
# The following functions  shift_axis_y()  and  shift_axis_x()
# are applied to place x and y axes in the plots
# to hazard ratios = "1" (instead of "0")

shift_axis_y <- function(p, y=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(y=y)
  ax <- g[["grobs"]][g$layout$name == "axis-b"][[1]]
  p + annotation_custom(grid::grobTree(ax, 
                                       vp = grid::viewport(
                                         y=0, height=sum(ax$height))), 
                        ymax=y, ymin=y) +
    geom_hline(aes(yintercept=y), data = dummy) +
    theme(axis.text.x = element_blank(), 
          axis.ticks.x=element_blank())
}

shift_axis_x <- function(p, x=0){
  g <- ggplotGrob(p)
  dummy <- data.frame(x=x)
  ax <- g[["grobs"]][g$layout$name == "axis-l"][[1]]
  p + annotation_custom(grid::grobTree(ax, 
                                       vp = grid::viewport(
                                         x=0, width = sum(ax$height))), 
                        xmax=x, xmin=x) +
    geom_vline(aes(xintercept=x), data = dummy) +
    theme(axis.text.y = element_blank(), 
          axis.ticks.y=element_blank(), 
          legend.position = "none")
}


#' * -----------------------------------------------------------------------


# IMPORT ####

# Normalized Enrichment Scores (NES) of pathways for iDFS prognosis
#  from biopsy and resect samples

df <-read.table("biopsy_resect_HR_basedGSEA_20092024.txt",header=TRUE, sep='\t')

# Define four quadrants Q1-Q4 for coloring
df$quadrant <- dplyr::case_when(df$Resect > 0 & df$Biopsy > 0  ~ "Q1", #Q1... both x and y are positive
                                df$Resect > 0 & df$Biopsy <0  ~ "Q2",  #Q2... x is neg and y is positive
                                df$Resect < 0 & df$Biopsy <0  ~ "Q3",  #Q3... both x and y are negative
                                df$Biopsy > 0 & df$Resect <0  ~ "Q4")  #Q4... x is pos and y is negative

testplot<-ggplot(df, aes(x = Biopsy, y =Resect, color = quadrant), guides(fill = FALSE, color = FALSE)) +
  geom_point(size=2) +
  geom_text_repel(aes(label = Pathway),size=3.5, box.padding = 0.5, max.overlaps=Inf) +
  theme(panel.background = element_rect(fill="white"),  
        plot.margin = margin(1, 1, 1, 1, "cm"),
        axis.title=element_blank(),
        axis.ticks.length=unit(.2, "cm"))+
  scale_x_continuous(limits = c(-3, 3),  n.breaks=7) +
  scale_y_continuous(limits = c(-3, 3),n.breaks=7)+
  labs(title = "paired 540 samples", 
       caption = "NES: normalized enrichment score \n(pos. NES: genes set enriched in tumors from patients with good prognosis; neg. NES = gene set enriched in tumors from patients with poor prognosis)") +
  theme(plot.title = element_text(hjust = 0, vjust = 10, margin = margin(t = 10)),
        plot.caption = element_text(hjust = 0))

# Add central axes and axis labels
  
p1 <- testplot +
  annotate("segment", x = -Inf, xend = Inf, y = 0, linewidth=1.5, yend = 0, color = "purple") +
  annotate("segment", x = 0, xend = 0, y = -Inf, linewidth=1.5, yend = Inf, color = "blue") +
  
  annotate("text", x = -Inf, y = 1, label = "NES for iDFS in post.Tx \u2192  improved survival", angle = 90,size=4, hjust=1.0, vjust=-2.0,  color = "red") +
  annotate("text", x = -1, y = Inf, label = "NES for iDFS in pre.Tx \u2192  improved survival", hjust=0,vjust=30, size=4, color = "blue") +
  
  annotate("text", x = -Inf, y = Inf, label = "good iDFS in post.Tx &\npoor iDFS in pre.Tx", vjust = 1, hjust=0, size=3.5)+
  annotate("text", x = -Inf, y = -Inf, label = "poor iDFS in post.Tx & pre.Tx", vjust = 0, hjust=0, size=3.5)+
  annotate("text", x = Inf, y = Inf, label = "good iDFS in post.Tx & pre.Tx", vjust = 1, hjust=1, size=3.5)+
  annotate("text", x = Inf, y = -Inf, label = "good iDFS in pre.Tx &\npoor iDFS in post.Tx", vjust = 0, hjust=1, size=3.5)+
  
  coord_cartesian(clip = "off")  # Allow annotations to be outside the plot area

# Shift axis labels  
 
p1<-shift_axis_y(p1, y=0)
p1<-shift_axis_x(p1, x=0)
p1

# Save plot to file
  
ggsave(file="SFig5.svg", plot=p1, width=10, height=8)
dev.off()


#'\newpage
# SESSION INFO ####
sessionInfo()
