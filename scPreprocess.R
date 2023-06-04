
scPreprocess = function(data_path, cell_call_dir = "",
                      experiment_group = "UW44scx01", project="Project"){

  ###################################################
  #USING SEURAT TO PROCESS RAW DATA
  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )
  if (!require("sctransform", quietly = TRUE))
    install.packages("sctransform", library(sctransform) )
  if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse", library(tidyverse) )
  if (!require("glmGamPoi", quietly = TRUE))
    install.packages("glmGamPoi", library(glmGamPoi) )
  if (!require("harmony", quietly = TRUE))
    install.packages("harmony", library(harmony) )

  if(cell_call_dir != ""){
    cell.calls<-load(cell_call_dir)
    experiment_call<-as.data.frame(Cells.class[experiment_group])
    colnames(experiment_call)<-NULL
  }
  else{
    #NO CELL CALL PROVIDED.
    #SKIPPING CELL CALL PREPROCESSING
    print()
  }
  # library(Seurat)
  # library(SeuratData)
  # library(BPCells)
  # library(dplyr)
  options(Seurat.object.assay.version = "v5")

  min.cells=10
  min.features=200
  sample.data<-Read10X(data_path)
  cell.names<-colnames(sample.data)
  cell.names<-sub("-1", "", cell.names)
  colnames(sample.data)<-cell.names

  seurat_obj<-CreateSeuratObject(sample.data,project=project,min.cells=min.cells, min.features=min.features)
  seurat_obj<-PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.MT")

  if(cell_call_dir != "")
    seurat_obj$cell.calls<-experiment_call
    table(seurat_obj$cell.calls)
    Idents(seurat_obj)<-'cell.calls'
    seurat_obj$NA.calls<-is.na(seurat_obj$cell.calls)
    seurat_obj<-subset(seurat_obj, NA.calls=='FALSE')

  plot1<-VlnPlot(seurat_obj,features=c('nCount_RNA','nFeature_RNA','percent.MT'))
  plot1

  return (seurat_obj)
}
