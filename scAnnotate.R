scAnnotate = function(object, cell_matrix = 'SCT', tissue = 'Brain'){
  if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", library(dplyr) )
  if (!require("openxlsx", quietly = TRUE))
    install.packages("openxlsx", library(openxlsx) )
  if (!require("HGNChelper", quietly = TRUE))
    install.packages("HGNChelper", library(HGNChelper) )

  ###Source and database
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"; tissue = c(tissue,"Immune system")
  gs_list = gene_sets_prepare(db_, tissue)

  #Set gs2 to NULL if non negative markers
  es.max = sctype_score(scRNAseqData = object[[cell_matrix]]@scale.data,
                        scaled = TRUE, gs = gs_list$gs_positive,
                        gs2 = gs_list$gs_negative)

  #Merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(object@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(object@meta.data[object@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(object@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)


  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/6.7] = "Unknown"
  #Sctype score
  print(sctype_scores[,1:3])

  #Plot
  object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    object@meta.data$customclassif[object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }


  plot <- DimPlot(object, reduction = "umap", label = TRUE, repel = TRUE,
                  group.by = 'customclassif')
  plot

  return (object)

}
