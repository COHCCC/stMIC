library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(Matrix)
library(e1071)
library(ggspatial)
library(viridis)

#Note: The file 'relabelled_clusters_R_E.csv' is the output from the tutorial. 
# It's a combined feature matrix (ResNet-gene-Expression) with a clustering label for each barcode. 
# 'Relabelled' indicates that it has been harmonized to maximize overlap with original Louvain cluster results 
# for improved visualization and comparison.

# working space
setwd('../sample')
# the following code is to merge relabelled combined method clusters to the original table.
file1  <- read.table("./annotation/relabelled_clusters_R_E.csv", sep = ",", header = TRUE)
file2 <- read.table("./spatial/tissue_positions_list.csv", sep = ",", header = FALSE)

# Perform partial matching of barcodes
matched_rows <- sapply(file1$barcode, function(x) any(grepl(x, file2$V1)))
matched_data <- file2[matched_rows, ]
merged_data <- merge(file1, matched_data, by.x = "barcode", by.y = "V1", all.y = TRUE)
merged_data <- merged_data[, c(colnames(merged_data)[colnames(merged_data) != "cluster"], "cluster")]
merged_data[is.na(merged_data)] <- ""
# View the merged data
print(merged_data)
write.csv(merged_data, "./annotation/relabelled_RE_Visulization.csv", quote = FALSE, row.names = FALSE, col.names = FALSE)
# remember to remove the first row of this output file for the following analysis


geom_spatial <-  function(mapping = NULL,
                          data = NULL,
                          stat = "identity",
                          position = "identity",
                          na.rm = FALSE,
                          show.legend = NA,
                          inherit.aes = FALSE,
                          ...) {
  
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


sample_names <- c("Combined clustering result (ResNet-gene-Expression)")

tissue_paths <- c("./annotation/relabelled_RE_Visulization.csv")
image_paths <- c("./spatial/tissue_hires_image.png") # been provided in "spatial folder"
scalefactor_paths <- c('./spatial/scalefactors_json.json') # been provided in spatial folder"

images_cl <- list()
for (i in 1:length(sample_names)) {
  images_cl[[i]] <- read.bitmap(image_paths[i])
}

height <- list()

for (i in 1:length(sample_names)) {
  height[[i]] <-  data.frame(height = nrow(images_cl[[i]]))
}

height <- bind_rows(height)

width <- list()

for (i in 1:length(sample_names)) {
  width[[i]] <- data.frame(width = ncol(images_cl[[i]]))
}

width <- bind_rows(width)


grobs <- list()
for (i in 1:length(sample_names)) {
  grobs[[i]] <- rasterGrob(images_cl[[i]], width=unit(1,"npc"), height=unit(1,"npc"))
}

images_tibble <- tibble(sample=factor(sample_names), grob=grobs)
images_tibble$height <- height$height
images_tibble$width <- width$width

scales <- list()

for (i in 1:length(sample_names)) {
  scales[[i]] <- rjson::fromJSON(file = scalefactor_paths[i])
}


bcs <- list()
for (i in 1:length(sample_names)) {
  bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol","Cluster"), header = FALSE)
  bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_hires_scalef    # scale tissue coordinates for highres image
  bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_hires_scalef
  bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
  bcs[[i]]$height <- height$height[i]
  bcs[[i]]$width <- width$width[i]
}

names(bcs) <- sample_names
bcs_merge <- bind_rows(bcs, .id = "sample")


## Plotting


plot_color=c("#a6cee3","#1f78b4","#b2df8a","#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00","#cab2d6","#997273","#787878","#db4c6c","#9e7a7a","#554236","#af5f3c","#93796c","#f9bd3f","#dab370","#877f6c","#268785")

plots <- list()

for (i in 1:length(sample_names)) {
  
  plots[[i]] <- bcs_merge %>%
    filter(sample ==sample_names[i]) %>%
    filter(tissue == "1") %>%
    ggplot(aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_manual(values = plot_color)+
    xlim(0,max(bcs_merge %>%
                 filter(sample ==sample_names[i]) %>%
                 select(width)))+
    ylim(max(bcs_merge %>%
               filter(sample ==sample_names[i]) %>%
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    labs(fill = "Cluster")+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    theme_set(theme_bw(base_size = 12))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.text = element_blank(),
          axis.ticks = element_blank())
}

plot_grid(plotlist = plots)

for (i in 1:length(sample_names)) {
  plots[[i]] <- bcs_merge %>%
    filter(sample ==sample_names[i]) %>%
    filter(tissue == "1") %>%
    ggplot(aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
    geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
    geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
    coord_cartesian(expand=FALSE)+
    scale_fill_manual(values = plot_120D)+
    xlim(0,max(bcs_merge %>%
                 filter(sample ==sample_names[i]) %>%
                 select(width)))+
    ylim(max(bcs_merge %>%
               filter(sample ==sample_names[i]) %>%
               select(height)),0)+
    xlab("") +
    ylab("") +
    ggtitle(sample_names[i])+
    labs(fill = "Cluster")+
    guides(fill = guide_legend(override.aes = list(size=4)))+
    theme_set(theme_bw(base_size = 12))+
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent"), # Makes plot background transparent
      axis.line = element_line(colour = "black"),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
}


##########################
# NES table
##########################
####################################
# Side by side bar plot
####################################
library(ggplot2)
data <- data.frame(Group = c("combined", "louvain"),
                   Count = c(202, 105))

ggplot(data, aes(x = Group, y = Count, fill = Group)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = c("#1f78b4", "#db4c6c")) +
  labs(x = "",y = "Count", title = "Number of Significant Pathways") +
  theme_minimal()

df <- read.csv("/Users/ninasong/Desktop/spatialProject/spatialDataset/slide117-120/slide120/slide120_D1/boxplot_4.csv")

library(tidyverse)

# Create a long format of the dataframe
df_long <- df %>%
  pivot_longer(everything(), names_to = "Group", values_to = "Value")
# Identify the top three outliers per group
df_outliers <- df_long %>%
  group_by(Group) %>%
  top_n(3, Value)
# Create a boxplot
p <- ggplot(df_long, aes(x = Group, y = Value, fill = Group, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.3, lwd = 1.5) +
  scale_fill_manual(values = c("lightblue", "#fb9a99")) +
  scale_color_manual(values = c("#1f78b4", "#db4c6c")) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = "none"
  ) +
  labs(x = "", y = "Normalized Enrichment Score")

p


#################################
# estimate analysis
#################################
library(ggplot2)
value = c("#a6cee3","#b2df8a", "#ff7f00","#33a02c","#1f78b4","gold", "#db4c6c", "#a65628","#4daf4a", "grey","purple")
setwd('/Users/ninasong/Desktop/spatialProject/spatialDataset/slide117-120/slide120/slide120_D1/estimateScore/')
stromalScore_ge <- read.csv(file = "stromalScore_RE.csv", sep = ",")
stromalScore_ge$Cluster <- as.factor(stromalScore_ge$Cluster)
# Create the plot
p <- ggplot(stromalScore_ge, aes(x = Cluster, y = Score, fill = Cluster)) + 
  geom_violin(trim=FALSE, scale="area", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  scale_fill_manual(values=value) +
  labs(title="Stormal score for combined method", x="Cluster", y="Score") +
  theme_classic() +
  theme(text = element_text(size=15), # increase base font size
        panel.grid.major = element_blank(), # remove major grid
        panel.grid.minor = element_blank(), # remove minor grid
        axis.line = element_line(colour = "black"), # add solid axes lines
        legend.position = "none", # remove legend
        plot.title = element_text(hjust = 0.5)) # rotate x-axis text

# Display the plot
print(p)


################################ 
#install CePa package for read gct functionality
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CePa")
library(CePa)
setwd("../sample")

## manually add cluster labels to Seurat Object metadata
ssgsea <- as.data.frame(read.gct("../FFD1_expression.gct"))

ST.orig <- Load10X_Spatial(data.dir = ".")
#get barcodes from spatial object and create vector, then replace names in ssgsea dataframe
barcodes <- c(Cells(ST.orig))
colnames(ssgsea) <- barcodes

#create new assay to store ssGSEA information in the Surat object
assay <- CreateAssayObject(counts = ssgsea)
ST.orig[["ssGSEA"]] <- assay
#remove(ssgsea, assay)

DefaultAssay(FFD1) <- "ssGSEA"
ST.orig <- SCTransform(ST.orig, assay = "ssGSEA", verbose = FALSE)

ST.meta <- read.csv("./annotation/relabelled_clusters_RE.csv", sep = ',')
ST.meta <- `colnames<-`(ST.meta, c("barcode", "cluster"))
rownames(ST.meta) <- ST.meta$barcode
## check metadata
head(ST.meta)

ST <- AddMetaData(object = ST.orig, metadata = ST.meta)
head(ST@meta.data)
Idents(ST = ST) <- "cluster"

plot_120D = c("#a6cee3","#b2df8a","#ff7f00","#33a02c","#1f78b4", "gold","#e31a1c","#787878","#ff7f00","#fb9a99","#1f78b4","#e31a1c")

SpatialPlot(ST, group.by = "cluster", crop = FALSE, image.alpha = 0,pt.size.factor = 1.32) +
  theme(panel.grid = element_blank()) +
  scale_fill_manual(values = plot_120D)

SpatialPlot(ST.orig, features = c("IFIT1"), crop = FALSE, image.alpha = 0, pt.size.factor = 1.33) +
  theme(panel.grid = element_blank()) +
  scale_fill_distiller(palette = "Spectral")
