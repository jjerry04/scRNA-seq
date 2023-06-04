
scBubble = function(object, cell_matrix = 'SCT', tissue = 'Brain',
                    reduction ='umap', labelSize = 1.5, score = 4,
                    title = "HVG"
                    ){

  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", library(BiocManager) )
  if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", library(dplyr) )
  if (!require("openxlsx", quietly = TRUE))
    install.packages("openxlsx", library(openxlsx) )
  if (!require("HGNChelper", quietly = TRUE))
    install.packages("HGNChelper", library(HGNChelper) )
  if (!require("igraph", quietly = TRUE))
    install.packages("igraph", library(igraph))
  if (!require("data.tree", quietly = TRUE))
    install.packages("data.tree", library(data.tree))
  if (!require("ggraph", quietly = TRUE))
    install.packages("ggraph", library(ggraph))
  if (!require("tidyverse", quietly = TRUE))
    install.packages("tidyverse", library(tidyverse))

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
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores*0.5)


  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/score] = "Unknown"

  #Plot
  object@meta.data$customclassif = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    object@meta.data$customclassif[object@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }

  # prepare edges
  cL_resutls=cL_resutls[order(cL_resutls$cluster),];
  edges = cL_resutls;
  edges$type = paste0(edges$type,"_",edges$cluster);
  edges$cluster = paste0("cluster ", edges$cluster);
  edges = edges[,c("cluster", "type")];
  colnames(edges) = c("from", "to");
  rownames(edges) <- NULL

  # prepare nodes
  nodes_lvl1 = sctype_scores[,c("cluster", "ncells")];
  nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster);
  nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1;
  nodes_lvl1$realname = nodes_lvl1$cluster;
  nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c();
  ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06",
                     "#eccf5a","#b5aa0f","#e4b680","#7ba39d",
                     "#b15928","#ffff99", "#6a3d9a","#cab2d6",
                     "#ff7f00","#fdbf6f","#e31a1c","#fb9a99",
                     "#33a02c","#b2df8a","#1f78b4","#a6cee3")
  for (i in 1:length(unique(cL_resutls$cluster))){
    dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ];
    nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster =
                                                paste0(dt_tmp$type,"_",dt_tmp$cluster),
                                              ncells = dt_tmp$scores, Colour = ccolss[i],
                                              ord = 2, realname = dt_tmp$type))
  }

  nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
  files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")];
  files_db = unique(files_db); nodes = merge(nodes, files_db,
                                             all.x = T, all.y = F,
                                             by.x = "realname",
                                             by.y = "cellName", sort = F)
  nodes$shortName[is.na(nodes$shortName)] =
    nodes$realname[is.na(nodes$shortName)];
  nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

  for( i in unique(nodes$realname)){
    a<-subset(nodes, realname == i)
    if( length(unique(a$shortName)) > 1){
      nodes[nodes["realname"] == i,"shortName"]<-min(char(unique(a$shortName)))
    }
  }
  nodes<-nodes[!duplicated(nodes),]

  mygraph <- graph_from_data_frame(edges, vertices=nodes)

  # Make the graph
  gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) +
    geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) +
    geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
    theme_void() +
    geom_node_text(aes(filter=ord==2, label=shortName,
                                      colour=I("black"), fill="white",
                                      repel = !1, parse = T,
                                      size = I(log(ncells,25)*labelSize)))
                                      # + geom_node_label(aes(filter=ord==1,
                                      # label=shortName, colour=I("#000000"),
                                      # size = I(2), fill="white",
                                      # parse = T), repel = !0, segment.linetype="dotted")


  print(head(object))
  plot <- DimPlot(object, reduction = reduction,
          label = TRUE,
          cols = ccolss, repel = TRUE)

  bubbleplot <- scater::multiplot(plot, gggr +ggtitle(title) , cols = 2)
  #bubbleplot <- scater::multiplot(bubbleplot, gggr, cols = 2)
  return (bubbleplot)
}


scAutoTissue <- function(object,assay ='SCT', scale = TRUE ){
  # load auto-detection function
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/auto_detect_tissue_type.R")
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

  # guess a tissue type
  tissue_guess = auto_detect_tissue_type(path_to_db_file =
                                           db_, seuratObject = object,
                                         scaled = scale, assay = assay)  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used
  return (tissue_guess)
  }
