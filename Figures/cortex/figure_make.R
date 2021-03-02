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

genespaceARIs = read.csv('genespaceARIs.csv')
genespaceARIs = rbind(genespaceARIs, read.csv('genespaceARIs_competitors.csv'))[-1]
raw = genespaceARIs[genespaceARIs$Method == "Raw",]
genespaceARIs = genespaceARIs[genespaceARIs$Method != "Raw",]
getlink = function(method, type) paste0(method, "_", type, '.csv')

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

UMAPs1 = list()

link = getlink(methods[1], type = 'HVG')
UMAP_data = read.csv(link, row.names = NULL)[-1]

nb.cols <- length(unique(UMAP_data$Cell.Type))
mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

tmp = genespaceARIs[genespaceARIs$Type=='HVG',]
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs1[[i]] = getplots(UMAP_data) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
}

plt = getplots(UMAP_data, make_legend = T) + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
  scale_fill_manual(values = mycolors3)
ggsave("type_key.tiff", plt, dpi = 300)

UMAPs2 = list()
for(i in c(1:length(methods))) {
  link = getlink(methods[i], type = 'HVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  UMAPs2[[i]] = getplots(UMAP_data, batch_indicator = 'Sequencing.Method')  + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
}

plt = getplots(UMAP_data, make_legend = T, batch_indicator = 'Sequencing.Method') + ggtitle(paste0(methods[i], " (HVGs, ARI=", round(tmp$ARI[i], digits = 2), ")")) +
  scale_color_brewer(palette = "Set2")
ggsave("type_batch.tiff", plt, dpi = 300)

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
  UMAPs4[[i]] = getplots(UMAP_data, batch_indicator = 'Sequencing.Method') + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
}

full = plot_grid(plotlist = c(UMAPs1, UMAPs2,UMAPs3, UMAPs4), ncol = sum(genespaceARIs$Type == 'HVG'))

ggsave('main.tiff', height = 5 * 1.1, width = 7 * 1.1, full, dpi = 300, units = "in")

ARIs = read.csv('../cortex_embedding/ARIsummary.csv')
ARIs = rbind(ARIs, read.csv('../cortex_embedding/ARIsummary_competitors.csv'))[-1]
ARIs = rbind(ARIs, raw[raw$Type == 'HVG',][-ncol(raw)])
ARIs <- gather(ARIs, Metric, Score, c('ARI', 'NMI', 'Purity')) #Create long format
ARIs$Metric = as.factor(ARIs$Metric)  
ARIs$Method = as.factor(ARIs$Method)

myplot <- ggplot(data = ARIs)
myplot <- myplot + geom_bar(aes(Method, Score, fill = Metric), stat="identity", position=position_dodge()) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
myplot = myplot + scale_fill_brewer(palette="Pastel2") + theme(text = element_text(size = 10)) + ggtitle('Clustering Accuracy with Embedding')
myplot = myplot + theme(plot.title = element_text(size=12)) + theme(plot.title = element_text(hjust=0.5))
myplot

ggsave('embeddingplot.tiff', height = 3.65, width = 6.35, units = "in", myplot, dpi = 300)

