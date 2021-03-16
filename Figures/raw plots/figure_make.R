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

datasets = c("Mouse Retina",
             "Macaque Retina",
             "Mouse Cortex",
             "Pancreas",
             "PBMC")

paths = c("../mouse_retina",
          "../macaque",
          "../cortex",
          "../pancreas_supplement",
          "../PBMC")

Batch_IDs = c("BatchID", "Macaque", "Sequencing.Method", "Technology", "Sequencing.Method")

getlink = function(path, type) paste0(path, '/Raw_', type, '.csv')
getARI = function(path, type) { 
  tabel = read.csv(paste0(path, '/Raw_ARIs.csv'))
  return(tabel$ARI[tabel$Type == type])
}
getplots = function(data, batch_indicator = NA){
  if(is.na(batch_indicator)) {
    myplot = 
      ggplot(data) + 
      geom_point(aes(x = UMAP1, y = UMAP2, color = Cell.Type), shape='.', show.legend = FALSE) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  } else{
    myplot =
      ggplot(data) + 
      geom_point(aes_string(x = "UMAP1", y = "UMAP2", color = batch_indicator), shape='.', show.legend = FALSE) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
  }
  
  myplot = myplot + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x =element_blank(),
                          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y =element_blank())
  
  myplot = myplot + theme(plot.title = element_text(hjust=0.5)) + theme(plot.title = element_text(size=18))
  
  return(myplot)
}

for(i in c(1:length(paths))) {
  dir.create(datasets[i])
  UMAPs = list()
  
  link = getlink(paths[i], type = 'HVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  
  nb.cols <- length(unique(UMAP_data$Cell.Type))
  mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)

  HVG_ARI = getARI(paths[i], "HVG")
  myplot = getplots(UMAP_data) + ggtitle(paste0(datasets[i], " (HVGs, ARI=", round(HVG_ARI, digits = 2), ")")) +
  scale_fill_manual(values = mycolors3)
  
  key = "Cell.Type"
  variable = unique(UMAP_data[,colnames(UMAP_data) == key])
  
  nulldf_ = data.frame(as.numeric(rep(NA, length(variable))))
  
  for(j in c(2:ncol(UMAP_data))){
    nulldf_[,j] = as.numeric(rep(NA, length(variable)))
  }
  
  colnames(nulldf_) = colnames(UMAP_data)
  nulldf_[,colnames(nulldf_) == key] = as.factor(variable)
  
  HVGcellplot = myplot
  
  HVGcellplot_legend = myplot + geom_point(data = nulldf_, aes(x = UMAP1, y = UMAP2, color = Cell.Type)) + 
    theme(legend.position = "right", legend.title = element_blank()) + theme(legend.key = element_rect(fill = "white"))
  
  key = Batch_IDs[i]
  myplot = getplots(UMAP_data, batch_indicator = key)  + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
  
  variable = unique(UMAP_data[,colnames(UMAP_data) == key])
  
  nulldf_ = data.frame(as.numeric(rep(NA, length(variable))))
  
  for(j in c(2:ncol(UMAP_data))){
    nulldf_[,j] = as.numeric(rep(NA, length(variable)))
  }
  
  colnames(nulldf_) = colnames(UMAP_data)
  nulldf_[,colnames(nulldf_) == key] = as.factor(variable)
  
  HVGbatchplot = myplot
  HVGbatchplot_legend = myplot + geom_point(data = nulldf_, aes_string(x = "UMAP1", y = "UMAP2", color = key)) + 
    theme(legend.position = "right", legend.title = element_blank()) + theme(legend.key = element_rect(fill = "white"))
  
  link = getlink(paths[i], type = 'LVG')
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  
  nb.cols <- length(unique(UMAP_data$Cell.Type))
  mycolors3 <- colorRampPalette(brewer.pal(12, "Set3"))(nb.cols)
  UMAP_data = read.csv(link, row.names = NULL)[-1]
  
  LVG_ARI = getARI(paths[i], "LVG")
  myplot = getplots(UMAP_data) + ggtitle(paste0(datasets[i], " (LVGs, ARI=", round(LVG_ARI, digits = 2), ")")) +
    scale_fill_manual(values = mycolors3)
  
  key = "Cell.Type"
  variable = unique(UMAP_data[,colnames(UMAP_data) == key])
  
  LVGcellplot = myplot
  
  key = Batch_IDs[i]
  myplot = getplots(UMAP_data, batch_indicator = key)  + ggtitle(" ") +
    scale_color_brewer(palette = "Set2")
  
  variable = unique(UMAP_data[,colnames(UMAP_data) == key])
  
  LVGbatchplot = myplot
  
  full = plot_grid(HVGcellplot, LVGcellplot,  HVGbatchplot, LVGbatchplot, ncol = 2)
  
  ggsave(paste0(datasets[i], "/type_key.tiff"), dpi = 300, HVGcellplot_legend)
  ggsave(paste0(datasets[i], "/type_batch.tiff"), dpi = 300, HVGbatchplot_legend)
  ggsave(paste0(datasets[i], "/main.tiff"), height = 7.4, width = 8.5, units = "in", dpi = 300, full)
}