# This Rscript profiles the runtime of MNN

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
library(profmem)

path = "../Data/liver"
mydata = read_cortex(path)

meta_data = mydata[[2]]
genes = mydata[[3]]
mydata = mydata[[1]]
rownames(mydata) = rownames(meta_data)
colnames(mydata) = genes

fracs = c(0.1,0.2,0.4,0.6,0.8,1.0)

run_mnn = function(ds, md) {
  sf = librarySizeFactors(ds)
  countmatrix = normalizeCounts(ds, size_factors = sf)
  countmatrix = Matrix::as.matrix(countmatrix)
  batch_vec = md$sampleid
  nbatch = length(unique(batch_vec))
  
  return(mnnCorrect(countmatrix, batch = as.factor(batch_vec), auto.order = c(1:nbatch)))
}

profiled_mnn = data.frame(c(0), c(0), c("MNN"), c(as.integer(0)))

for(frac in fracs) {
  indices = read.csv(paste0("indices", frac, ".csv"))$X0
  subdata = mydata[, indices]
  meta_datatmp = meta_data[indices,]
  nCells = apply(subdata, MARGIN = 1, function(x) sum(x > 0) >= 1)
  subdata = subdata[nCells, ]
  start = Sys.time()
  print(start)
  mem_used = total(profmem({run_mnn(subdata, meta_datatmp)}))/1048576
  time_used = as.numeric(Sys.time() - start, units = "secs")
  profiled_mnn[nrow(profiled_mnn) + 1,] = c(time_used, mem_used, "MNN", as.integer(100*frac))
  tmp = profiled_mnn
  colnames(tmp) = c("Time (Seconds)", "Memory (MiB)", "Method", "Percent")
  write.csv(profiled_mnn[-1,], "../Figures/liver/mnn_profile.csv")
}


colnames(profiled_mnn) = c("Time (Seconds)", "Memory (MiB)", "Method", "Percent")
profiled_mnn