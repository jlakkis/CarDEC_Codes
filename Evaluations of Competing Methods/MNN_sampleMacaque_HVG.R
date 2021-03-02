setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

if("Seurat" %in% rownames(installed.packages()) == F) {
  install.packages('Seurat')
}

if("reticulate" %in% rownames(installed.packages()) == F) {
  install.packages('reticulate')
}

if("batchelor" %in% rownames(installed.packages()) == F) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("batchelor")
}

if("scater" %in% rownames(installed.packages()) == F) {
  BiocManager::install("scater")
}

reticulate::use_condaenv("cardec", required = T)
reticulate::source_python("mnn_reticulate.py")

library(Seurat)
library(batchelor)
library(scater)
library(Matrix)

path = "../Data/macaque_bc.h5ad"
mydata = read_macaque(path)

meta_data = mydata[[2]]
genes = mydata[[3]]
mydata = mydata[[1]]
rownames(mydata) = rownames(meta_data)
colnames(mydata) = genes

mydata = CreateSeuratObject(t(mydata))

hvgs = get_hvgs(mydata@assays$RNA@data, meta_data$sample, rownames(mydata@assays$RNA@data))
mydata = subset(mydata, features = hvgs)

sf = librarySizeFactors(mydata@assays$RNA@counts)
countmatrix = normalizeCounts(mydata@assays$RNA@counts, size_factors = sf)
countmatrix = Matrix::as.matrix(countmatrix)
batch_vec = meta_data$sample
nbatch = length(unique(batch_vec))

start = Sys.time()

output = mnnCorrect(countmatrix, batch = meta_data$sample, merge.order = c(1:nbatch))

Sys.time() - start

path = "MNNcorrected_hvg"
dir.create(path)

write.csv(output@assays@data$corrected, paste0(path, "/corrected_dataSample.csv"))
write.csv(meta_data, paste0(path, "/corrected_metadataSample.csv"))