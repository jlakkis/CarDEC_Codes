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

path = "../Data/cortex"
mydata = read_cortex(path)

meta_data = mydata[[2]]
genes = mydata[[3]]
mydata = mydata[[1]]
rownames(mydata) = rownames(meta_data)
colnames(mydata) = genes

mydata = CreateSeuratObject(t(mydata))

hvgs = get_hvgs(mydata@assays$RNA@data, meta_data$Method, rownames(mydata@assays$RNA@data))
mydata = subset(mydata, features = hvgs)

sf = librarySizeFactors(mydata@assays$RNA@counts)
countmatrix = normalizeCounts(mydata@assays$RNA@counts, size_factors = sf)
countmatrix = Matrix::as.matrix(countmatrix)
batch_vec = meta_data$Method
nbatch = length(unique(batch_vec))

start = Sys.time()

output = mnnCorrect(countmatrix, batch = as.factor(batch_vec), merge.order = c(1:nbatch))

Sys.time() - start

path = "MNNcorrected_hvg"
dir.create(path)

write.csv(output@assays@data$corrected, paste0(path, "/corrected_data_Cortex.csv"))
write.csv(meta_data, paste0(path, "/corrected_metadata_Cortex.csv"))