setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(png))
suppressPackageStartupMessages(library(RANN))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(kableExtra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(vioplot))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(tidyr))

UMAP_data = read.csv('CVs_nclust.csv')
getlink = function(method, type) paste0(method, '.csv')

clust_nums = (ncol(UMAP_data) - 3)/3
plots = as.list(c(1:clust_nums))

for(i in c(1:clust_nums)) {
  index = 2 + 2 * i
  data = UMAP_data[,c(3, index, index + 1)]
  colnames(data) = c("Cell.Type", "UMAP1", "UMAP2")
  
  x = colnames(UMAP_data)[index]
  nclust = substr(x, start = gregexpr("_", x)[[1]][1] + 1, stop = nchar(x))
  
  myplot = 
    ggplot(data) + 
    geom_point(aes(x = UMAP1, y = UMAP2, color = Cell.Type), shape='.') + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  myplot = myplot + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x =element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y =element_blank())
  
  myplot = myplot + theme(legend.position = "none") + theme(plot.title = element_text(hjust=0.5)) + theme(plot.title = element_text(size=18))
  
  plots[[i]] = myplot + ggtitle(paste0(nclust, " Clusters by Type"))
  
  if(i==1){
  myplot$layers[[1]]$aes_params$shape=19
  myplot = myplot + theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.key = element_rect(fill = "white"))
  ggsave("type_key.tiff", myplot, dpi = 300)
  myplot$layers[[1]]$aes_params$shape="."
  }
}

plots2 = as.list(c(1:clust_nums))

for(i in c(1:clust_nums)) {
  index = 2 + 2 * i
  cluster_index = 3 + 2 * clust_nums + i
  data = UMAP_data[,c(cluster_index, index, index + 1)]
  colnames(data) = c("Clust_ID", "UMAP1", "UMAP2")
  data[,1] = sapply(data[,1], toString)
  
  x = colnames(UMAP_data)[index]
  nclust = substr(x, start = gregexpr("_", x)[[1]][1] + 1, stop = nchar(x))
  
  nb.cols <- as.numeric(nclust)
  mycolors <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
  
  myplot = 
    ggplot(data) + 
    geom_point(aes(x = UMAP1, y = UMAP2, color = Clust_ID), shape='.') + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  myplot = myplot + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x =element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y =element_blank())
  
  myplot = myplot + theme(legend.position = "none") + theme(plot.title = element_text(hjust=0.5)) + theme(plot.title = element_text(size=18))
  
  plots2[[i]] = myplot + ggtitle(paste0(nclust, " Clusters by Cluster")) +
    scale_fill_manual(values = mycolors)
  
  
  myplot$layers[[1]]$aes_params$shape = 19
  myplot = myplot + theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=5))) +
    theme(legend.key = element_rect(fill = "white"))
  ggsave(paste0("type_clust", nclust, ".tiff"), myplot, dpi = 300)
  myplot$layers[[1]]$aes_params$shape="."
}

full = plot_grid(plotlist = c(plots, plots2), ncol = 2, byrow = F)

ggsave('CVs_nclust.tiff', height = 10, width = 5, units = 'in', full, dpi = 300)