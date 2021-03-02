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

genespaceARIs = read.csv('../macaque/genespaceARIs.csv')
genespaceARIs = rbind(genespaceARIs, read.csv('../macaque/genespaceARIs_competitors.csv'))[-1]
raw = genespaceARIs[genespaceARIs$Method == "Raw",]
genespaceARIs = genespaceARIs[genespaceARIs$Method != "Raw",]
getlink = function(method, type) paste0('../macaque/', method, "_", type, '.csv')

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

methods = unique(genespaceARIs$Method)


link = getlink(methods[1], type = 'HVG')
UMAP_data = read.csv(link, row.names = NULL)[-1]

nb.cols <- length(unique(UMAP_data$Cell.Type))
mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

UMAPs1 = list()
tmp = genespaceARIs[genespaceARIs$Type=='HVG',]
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs1[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

plt = getplots(UMAP_data, make_legend = T) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")"))+
  scale_fill_manual(values = mycolors3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_key_region.tiff", plt, dpi = 300)

UMAPs2 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs2[[i]] = getplots(UMAP_data, batch_indicator = 'Region')  + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
}

plt = getplots(UMAP_data, make_legend = T, batch_indicator = 'Region') + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
  scale_color_brewer(palette = "Set2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_batch_region.tiff", plt, dpi = 300)

tmp = genespaceARIs[genespaceARIs$Type=='LVG',]

UMAPs3 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'LVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs3[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (LVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

UMAPs4 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'LVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs4[[i]] = getplots(UMAP_data, batch_indicator = 'Region') + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
}

full = plot_grid(plotlist = c(UMAPs1, UMAPs2,UMAPs3, UMAPs4), ncol = sum(genespaceARIs$Type == 'HVG'))

ggsave('main_region.tiff', height = 7.4, width = 8.5 * 5/4, units = 'in', full, dpi = 300)

ARIs = read.csv('ARIsummary.csv')
ARIs = rbind(ARIs, read.csv('ARIsummary_competitors.csv'))[-1]
ARIs = rbind(ARIs, raw[raw$Type == 'HVG',][-ncol(raw)])
ARIs <- gather(ARIs, Metric, Score, c('ARI', 'NMI', 'Purity')) #Create long format
ARIs$Metric = as.factor(ARIs$Metric)
ARIs$Method = as.factor(ARIs$Method)

############################################

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

methods = unique(genespaceARIs$Method)


link = getlink(methods[1], type = 'HVG')
UMAP_data = read.csv(link, row.names = NULL)[-1]

nb.cols <- length(unique(UMAP_data$Cell.Type))
mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
nb.cols <- length(unique(UMAP_data$Sample))
mycolors2 <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)

UMAPs1 = list()
tmp = genespaceARIs[genespaceARIs$Type=='HVG',]
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs1[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

plt = getplots(UMAP_data, make_legend = T) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")"))+
  scale_fill_manual(values = mycolors3) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_key_sample.tiff", plt, dpi = 300)

UMAPs2 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs2[[i]] = getplots(UMAP_data, batch_indicator = 'Sample')  + ggtitle(" ") +
    scale_fill_manual(values = mycolors2)
}

plt = getplots(UMAP_data, make_legend = T, batch_indicator = 'Sample') + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
  scale_fill_manual(values = mycolors2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
ggsave("type_batch_sample.tiff", plt, dpi = 300)

tmp = genespaceARIs[genespaceARIs$Type=='LVG',]

UMAPs3 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'LVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs3[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (LVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

UMAPs4 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'LVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs4[[i]] = getplots(UMAP_data, batch_indicator = 'Sample') + ggtitle(" ") +
    scale_fill_manual(values = mycolors2)
}

full = plot_grid(plotlist = c(UMAPs1, UMAPs2,UMAPs3, UMAPs4), ncol = sum(genespaceARIs$Type == 'HVG'))

ggsave('main_Sample.tiff', height = 7.4, width = 8.5 * 5/4, units = 'in', full, dpi = 300)

