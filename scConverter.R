#Convert from annData to seurat and seurat to
scConverter = function(data, datatype = "h5ad",
                       filename = ""){

  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )
  if (!require("SeuratData", quietly = TRUE))
    install.packages("SeuratData", library(SeuratData) )
  if (!require("SeuratDisk", quietly = TRUE))
    install.packages("SeuratDisk", library(SeuratDisk) )

  if (datatype == "h5ad" && filename != ""){
    SaveH5Seurat(data, filename = filename)
    filename = paste(filename, ".h5seurat")
    filename = str_replace_all(filename, pattern = " ", repl = "")
    Convert(filename, dest = datatype, verbose = TRUE, overwrite = TRUE)
  }else{
    print("You must sepcificy filename, dest, and datatype h5ad data type.")
    return (NULL)
  }

  if(datatype == "single_cell"){
    return (as.SingleCellExperiment(data))
  }
}



















