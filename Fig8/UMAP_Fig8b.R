# HEADER ####
#
# Version: 2024-11-14
#
# Figure 8B: UMAP color by AC-subtypes
# 
#
#
#
# SETUP ####
#'  
Sys.setenv(lang = "en_US")


#' *Install required packages if missing* -----------------------------------


# Package names from CRAN
packs <- c("ggplot2", "dplyr", "tibble", "ggnewscale", "umap")

# Install packages not yet installed
installed_packages <- packs %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packs[!installed_packages])
}

#' *Load required packages* -----------------------------------

invisible(lapply(packs, library, character.only = TRUE))


#'\newpage
# IMPORT ####

heatmap_br <- read.delim("UMAP335_Biop_Res.txt", stringsAsFactors = FALSE)
sampleinfo <- read.delim("UMAP_sample_info.txt", stringsAsFactors = TRUE)

heatmap_br$id <- NULL

# Analysis / Plots ####

filtered_expression_df <- t(heatmap_br)

set.seed(123)

umap_results <- umap(filtered_expression_df, n_neighbors = 15, min_dist = 0.3, metric = "euclidean")

kmeans_result <- kmeans(umap_results$layout, centers = 3)

umap_plot_df <- data.frame(umap_results$layout) %>%
  tibble::rownames_to_column("SampleName") %>%
  dplyr::inner_join(sampleinfo, by = "SampleName")%>%
  dplyr::mutate(Cluster = as.factor(kmeans_result$cluster) ) 



df_Bclusters <- data.frame(
  UMAP1 = umap_plot_df$X1,
  UMAP2 = umap_plot_df$X2,
  A.Cluster = umap_plot_df$AC1_5.clusters,
  pre.AIMS = umap_plot_df$pre.AIMS,
  AIMS = umap_plot_df$AIMS,
  Tissue=umap_plot_df$tissue,
  cls=umap_plot_df$Cluster
)

AC.colors <- c('AC.1' = "gold", 'AC.2' = "blue", 'AC.3' = 'red', 'AC.4' = "brown", 'AC.5' = "magenta")

# Create the plot
AC_ggplot <- ggplot(
  df_Bclusters,
  aes(
    x = UMAP1,
    y = UMAP2)
  ) + 
  geom_point(aes(color = A.Cluster, shape = Tissue), size=3) +
  scale_color_manual(name = "A.Cluster", values = AC.colors) +  
  new_scale_color() +  
  scale_shape_manual(name = "Tissue", values = c(17,19)) +  
  labs(color = "Legend") +  
  theme_bw() + 
  theme(
    axis.text = element_text(size = 12),  
    axis.title = element_text(size = 14),  
    legend.text = element_text(size = 12, family = "Arial"),
    legend.title = element_text(size = 12, family = "Arial")
  )

AC_ggplot
ggsave("UMAP_AC_clusters_Fig8B.svg", plot = AC_ggplot, device = "svg", width = 10, height = 8)

# SESSION INFO ####
sessionInfo()
