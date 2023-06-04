
# Pass objects to integrate
# exHarmony - exlcudes harmony batch correction
scIntegrate <- function(obj1, obj2, exHarmony = FALSE, 
  pattern = c('percent.MT','nCount_RNA')){

  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )
  if (!require("harmony", quietly = TRUE))
    install.packages("harmony", library(harmony) )

  iObject <- FindIntegrationAnchors(object.list =
                                      list(obj1, obj2),
                                    dims = 1:50
  )
  iObject <- IntegrateData(anchorset = iObject, dims = 1:50)
  head(iObject)

  #Run full downstream analysis
  iObject <- scDSA(object =iObject,
                       col_data = pattern,
                       reduct = 'pca', max_dim = 50, 
                       integrate = exHarmony)

  return (iObject)
}
