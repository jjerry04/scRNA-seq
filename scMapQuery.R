scMapQuery = function(reference, objectToQuery,
                      normalization.method = 'SCT',
                      reduction = 'pca',
                      dims = 50,
                      refData = list(),
                      reduction.model = 'umap'){

  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )

  anchors <- FindTransferAnchors( reference = reference, query = objectToQuery,
                                  normalization.method = normalization.method,
                                  reference.reduction = reduction, dims = 1:dims)

  #Construct REFERENCE MAP Querying
  mappedData <- MapQuery( anchorset = anchors, query = objectToQuery,
                     reference = reference, refdata = refData,
                     reference.reduction = reduction,
                     reduction.model = reduction.model)
  return (mappedData)
}
