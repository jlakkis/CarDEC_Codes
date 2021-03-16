setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()

packages = c("knitr", "cowplot", "gridExtra", "png", "RANN", "ggplot2", "kableExtra", "dplyr", "vioplot", "RColorBrewer", "tidyr")

for (package in packages) {
  if(package %in% rownames(installed.packages()) == F) {
    install.packages(package)
  }
}

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

ARIs = read.csv('ARIsummary.csv')
ARIs = rbind(ARIs, read.csv('ARIsummary_competitors.csv'))[-1]
methods = unique(ARIs$Method)

getlink = function(method, type) paste0(method, '.csv')

getplots = function(data, batch_indicator = NA, make_legend = F){
  if(is.na(batch_indicator)) {
    myplot = 
      ggplot(data) + 
      geom_point(aes(x = UMAP1, y = UMAP2, color = Cell.Type), shape='.') + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else{
    myplot =
      ggplot(data) + 
      geom_point(aes_string(x = "UMAP1", y = "UMAP2", color = batch_indicator), shape='.') + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  
  myplot = myplot + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x =element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y =element_blank())
  
  myplot = myplot + theme(legend.position = "none") + theme(plot.title = element_text(hjust=0.5)) + theme(plot.title = element_text(size=8))
  if(make_legend) {
    myplot = myplot + theme(legend.position = "right", legend.title = element_blank()) + guides(colour = guide_legend(override.aes = list(size=5))) +
      theme(legend.key = element_rect(fill = "white"))
    myplot$layers[[1]]$aes_params$shape=19
  }
  
  return(myplot)
}

link = getlink(methods[1], type = 'HVG')
UMAP_data = read.csv(link, row.names = NULL)[-1]

nb.cols <- length(unique(UMAP_data$Cell.Type))
mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

UMAPs1 = list()
tmp = ARIs
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs1[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

plt = getplots(UMAP_data, make_legend = T) + ggtitle(paste0(methods[i], " (ARI=", round(tmp$ARI[i], digits = 2), ")"))+
  scale_fill_manual(values = mycolors3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_key.tiff", plt, dpi = 300)

UMAPs2 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs2[[i]] = getplots(UMAP_data, batch_indicator = 'BatchID')  + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
}

plt = getplots(UMAP_data, make_legend = T, batch_indicator = 'BatchID') + ggtitle(paste0(methods[i], " (ARI=", round(tmp$ARI[i], digits = 2), ")")) +
  scale_color_brewer(palette = "Set2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_batch.tiff", plt, dpi = 300)

full = plot_grid(plotlist = c(UMAPs1, UMAPs2), ncol = length(UMAPs1))

ggsave('main.tiff', height = 7.4/2, width = 8.5 * 5/4, units = 'in', full, dpi = 300)

