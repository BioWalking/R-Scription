library(Seurat)
library(SeuratDisk)
library(Seurat)
fs = list.files('./GSE142425_RAW/',pattern = '^GSM')
samples <- substr(fs,1,10)
lapply( unique(samples), function(x){
  y = fs[grepl(x,fs)]
  folder = paste0('./GSE142425_RAW/',strsplit(y[1],split = '_')[[1]][1])
  #创建文件夹
  dir.create(folder,recursive = T)
  #重命名子文件夹并移动到相应的文件夹中
  file.rename(paste0('./GSE142425_RAW/',y[1]),file.path(folder,"barcodes.tsv.gz"))
  file.rename(paste0('./GSE142425_RAW/',y[2]),file.path(folder,"features.tsv.gz"))
  file.rename(paste0('./GSE142425_RAW/',y[3]),file.path(folder,"matrix.mtx.gz"))
})

samples <- dir(path="./GSE142425_RAW/", pattern="^GSM")
for (i in samples){
  assign(paste0( i), Read10X(data.dir = paste0("./GSE142425_RAW/", i)))
}

setwd("GSE142425_RAW")
sceList = lapply(samples,function(folder){ 
  CreateSeuratObject(min.cells = 5, min.features = 2000,counts = Read10X(folder), 
                     project = folder )##阈值
})

ss=c("E11","E13","E15","E18")
GSE142425 <- merge(sceList[[1]],
                 y = c(sceList[[2]],sceList[[3]],sceList[[4]]),
                 add.cell.ids = ss)
GSE142425$database="GSE142425"


setwd("..")
##GSE151985

Convert("GSE151985_annData_GLE_cells.h5ad", dest="h5seurat",
        assay = "RNA",
        overwrite=F)
GSE151985 <- LoadH5Seurat("GSE151985_annData_GLE_cells.h5seurat")
GSE151985 = UpdateSeuratObject(GSE151985)
GSE151985$database="GSE151985"
colnames(GSE151985@meta.data)[c(3:4,18)] <- c("nCount_RNA","nFeature_RNA","orig.ident")


##GSE179701
Convert("GSM5429626_ad_e17_raw_allcells.h5ad", dest="h5seurat",
        assay = "RNA",
        overwrite=F)
GSE179701 <- LoadH5Seurat("GSM5429626_ad_e17_raw_allcells.h5seurat")
GSE179701 = UpdateSeuratObject(GSE179701)
GSE179701$database="GSE179701"
GSE179701$orig.ident="e17.5"

Seurat_Data=merge(x=GSE142425,y=c(GSE151985,GSE179701))
save(Seurat_Data,file="Seurat_Data.Rdata")