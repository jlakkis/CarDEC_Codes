options(warn=-1) # turn off warning message globally
.libPaths(c("/home/xiaoxiang/R/x86_64-pc-linux-gnu-library/3.5",.libPaths()))
Sys.setenv(RETICULATE_PYTHON_ENV="/home/xiaoxiang/anaconda3/envs/py36")#="/home/xiaoxiang/.conda/envs/DESCVIR"
Sys.setenv(RETICULATE_PYTHON="/usr/bin/python3")
#RETICULATE_PYTHON="/home/xiaoxiang/anaconda3/bin/python3",
if ("Seurat" %in% loadedNamespaces()) detach("package:Seurat",unload = T)
dyn.load("/home/xiaoxiang/R/x86_64-pc-linux-gnu-library/3.5/sf/libs/sf.so")
#suppressPackageStartupMessages(library(monocle,lib.loc = "/usr/lib/R/monocle_alpha"))# devtools::install_github("")
#devtools::install_github("cole-trapnell-lab/DDRTree", ref="simple-ppt-like",lib="/usr/lib/R/monocle_alpha")
#devtools::install_github("r-spatial/sf") if 
#install.packages("~/Downloads/monocle-release-monocle3_alpha/", repos = NULL,lib = "/usr/lib/R/monocle_alpha")
suppressPackageStartupMessages(library(reticulate))
#suppressPackageStartupMessages(library(devtools))

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggjoy))
suppressPackageStartupMessages(library(VGAM))
suppressPackageStartupMessages(library(knitr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(kableExtra))
library(batchelor)
library(scater)
#library(Matrix)


Convert_to_seurat3=function(adata){
  suppressPackageStartupMessages(library("Seurat",lib.loc = "/usr/lib/R/self_library/"))
  mtx=py_to_r(adata$X$T$tocsc())
  cellinfo=py_to_r(adata$obs)
  geneinfo=py_to_r(adata$var)
  colnames(mtx)=cellinfo$cellname
  rownames(mtx)=rownames(geneinfo)
  obj=CreateSeuratObject(mtx,meta.data = cellinfo[,!colnames(cellinfo)%in%c("n_genes","n_counts"),drop=F],min.features  = 1)
  return(obj)
}

setwd("/home/xiaoxiang/Documents/carDEC_paper/CarDEC_new20200421/")
ad=import("anndata",convert = FALSE)
adata=ad$read_h5ad("../dca_test.h5ad")
obj0=Convert_to_seurat3(adata)
obj0=NormalizeData(obj0,verbose = F)
raw.data=obj0@assays$RNA@counts

sf = librarySizeFactors(raw.data)
countmatrix = normalizeCounts(raw.data, size_factors = sf)
# read highly variable genes CarDEC used
#df0=read.table("./CarDEC_hvg_used.tsv",sep="\t",header = T)
#gene_used=df0$genename[df0$highly_variable=="True"]
#countmatrix = Matrix::as.matrix(countmatrix[gene_used,])

countmatrix = Matrix::as.matrix(countmatrix)
batch_vec = obj0$dataset_batch
nbatch = length(unique(batch_vec))

start = Sys.time()

output = mnnCorrect(countmatrix, batch = as.factor(batch_vec), auto.order = c(1:nbatch))
colData(output)=DataFrame(cbind(as.data.frame(colData(output)),obj0@meta.data))
saveRDS(output,file="MNN_corrected_all.rds") #singlecell output@assays$data$corrected
#saveRDS(output,file="MNN_corrected.rds") #singlecell output@assays$data$corrected
Sys.time() - start
#path = "MNNcorrected"
#dir.create(path)
  