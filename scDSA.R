# Run downstream analysis on data. Set integrate to TRUE
# to create object for map query.
#single cell Downstream analysis (scDSA)
scDSA = function(object, col_data = c(), reduction = 'pca', assay = 'SCT',
                     group.by = 'orig.ident',max_dim = 50,
                     regress = c('percent.MT','nCount_RNA'), integrate = FALSE){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )
  if (!require("sctransform", quietly = TRUE))
    install.packages("sctransform", library(sctransform) )
  if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse", library(tidyverse) )
  if (!require("glmGamPoi", quietly = TRUE))
    BiocManager::install("glmGamPoi")
    library(glmGamPoi)
  if (!require("harmony", quietly = TRUE))
    install.packages("harmony", library(harmony) )

  #post integration
  if (integrate == FALSE){
    seurat.singlet<-SCTransform(object,vars.to.regress=regress,vst.flavor='v2')
    seurat.singlet<-RunPCA(seurat.singlet)
    head(seurat.singlet)
    seurat.singlet<-RunHarmony(seurat.singlet,assay.use=assay,
                               reduction=reduction,dims.use=1:max_dim,
                               group.by.vars=group.by,
                               plot.convergence=T)
    seurat.singlet<-FindNeighbors(seurat.singlet,reduction="harmony")
    seurat.singlet<-FindClusters(seurat.singlet,resolution=0.4)
    seurat.singlet<-RunUMAP(seurat.singlet,dims=1:50,reduction='harmony')
    plot1<-VlnPlot(seurat.singlet,features=c('nCount_RNA','nFeature_RNA','percent.MT'))
    plot1
    plot2 <- plot(DimPlot(seurat.singlet, label = TRUE))
    plot2
  }
  else{
    #Use for data integration (Do not cluster)
    seurat.singlet<-SCTransform(object,vars.to.regress=regress,vst.flavor='v2')
    head(seurat.singlet)
    seurat.singlet <- RunPCA(seurat.singlet)

    seurat.singlet<-RunHarmony(seurat.singlet,assay.use=assay,
                                reduction=reduction,dims.use=1:max_dim,
                                group.by.vars=group.by,
                                plot.convergence=T)
    seurat.singlet <- FindVariableFeatures(seurat.singlet,
                                    selection.method = "vst", nfeatures = 2000)
    seurat.singlet  <- ScaleData(seurat.singlet, assay = 'SCT')
    seurat.singlet <-RunUMAP(seurat.singlet,dims=1:50,reduction='pca',
                             return.model = TRUE )
    }

  return (seurat.singlet)
}
