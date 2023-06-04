
# -Log10 of genes(nFeature_RNA) detected/cell
# -log10 of UMI(nCounts_RNA)
# -divide by log10Genes/log10UMIs
# -Expect good cells above 0.8
# If there are many captured transcripts (high nUMI)
# and a low number of genes detected in a cell,
# this likely means that you only captured a low
# number of genes

nScore = function(object){
  #Return the novelty score of object
  #Determine RNA complexity through ratio of GENES over UMI
  genes = object$nFeature_RNA
  umis = object$nCount_RNA
  return(log10(genes)/log10(umis))
}


#   Calcuate proportion of transcript mapping
#   to MT genes. Instead of return percent value,
#   return the ratio value for pattern analysis.
#   Poor quality cells surpass 0.2 mt Ratio.

mitoRatio = function(object, pattern = '^MT-'){
  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )

  #Inverse percent values to obtain ratio values
  object$mtRatio = PercentageFeatureSet(object = object, pattern = pattern)
  object$mtRatio = object@meta.data$mtRatio/100
}

scQCVisual = function(object,include.all = TRUE, pattern = 'MT-',
                      countsPerSample = FALSE,
                      umiPerCell = FALSE,
                      genesPerCell = FALSE,
                      genesPerUMI = FALSE,
                      mtPerCell = FALSE,
                      quality = FALSE){

  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )
  if (!require("dplyr", quietly = TRUE))
    install.packages("dplyr", library(dplyr) )
  if (!require("magrittr", quietly = TRUE))
    install.packages('magrittr', library(magrittr))
  if (!require("ggplot2", quietly = TRUE))
    install.packages('ggplot2', library(ggplot2))

  metaData = object@meta.data
  cell_calls = object$cell.calls
  counts = object$nCount_RNA
  genes = object$nFeature_RNA

  #calculate nscore
  nscore = nScore(object)

  #Calculate pattern ratio
  patternRatio = mitoRatio(object, 'MT-')

  if(countsPerSample == TRUE){
    #Counts per sample
    if (include.all == TRUE){
      data <- metaData %>% ggplot(aes(x=sample, fill=cell_calls)) +
        geom_bar() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("Counts/Sample")
    }else{
      data <- metaData %>% ggplot(aes(x=msample, fill=cell_calls)) +
        geom_bar() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold")) +
        ggtitle("Counts/Sample")
    }

    return(data)
  }

  if(umiPerCell == TRUE){
    #Visualize UMIs
    # Visualize the number UMIs/transcripts per cell
    if (include.all == TRUE){
      data <- metaData %>%
        ggplot(aes(color=sample, x=counts, fill= sample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('UMI/Cell')
    }else{
      data <- metaData %>%
        ggplot(aes(color=msample, x=counts, fill= msample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        ylab("Cell density") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('UMI/Cell')
    }

    return (data)
  }
  if (genesPerCell == TRUE){
    # Visualize the distribution of genes detected per cell via histogram
    if(include.all == TRUE){
      data <- metaData %>%
        ggplot(aes(color=sample, x=genes, fill= sample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        scale_x_log10() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('Genes/Cell')
    }else{
      data <- metaData %>%
        ggplot(aes(color=msample, x=genes, fill= msample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        scale_x_log10() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('Genes/Cell')

    }
    return (data)
  }
  if (genesPerUMI == TRUE){
    # Visualize the overall complexity of the gene
    #expression by visualizing the genes detected per UMI (novelty score)
    if (include.all == TRUE){
      genesPerUMI = nScore(object)
      data <- metaData %>%
        ggplot(aes(x=genesPerUMI, color = sample, fill=sample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('Gene Complexity')
    }else{
      genesPerUMI = nScore(object)
      data <- metaData %>%
        ggplot(aes(x=genesPerUMI, color = msample, fill=msample)) +
        geom_density(alpha = 0.2) +
        theme_classic() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('Gene Complexity')
    }
    return (data)
  }
  if (mtPerCell == TRUE){
    # Visualize the distribution of mitochondrial gene expression detected per cell
    #Is there a large amount of MT contamination from dead
    #cells or dying cells
    if(include.all == TRUE){
      mtratio = mitoRatio(object)
      data <- metaData %>%
        ggplot(aes(color=sample, x=mtratio, fill=sample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('MT detected per cell')
    }else{
      mtratio = mitoRatio(object)
      data <- metaData %>%
        ggplot(aes(color=msample, x=mtratio, fill=msample)) +
        geom_density(alpha = 0.2) +
        scale_x_log10() +
        theme_classic() +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))+
        ggtitle('MT detected per cell')
    }
    return (data)

  }
  if(quality == TRUE){
    # Visualize the correlation between genes
    #detected and number of UMIs and determine
    #whether strong presence of cells with low
    #numbers of genes/UMIs exist.
    #Likely subset data in btm left quadrant
    #as

    mtratio = mitoRatio(object)
    data <- metaData %>%
      ggplot(aes(x=counts, y=genes, color=mtratio)) +
      geom_point() +
      scale_colour_gradient(low = "gray90", high = "black") +
      stat_smooth(method=lm) +
      scale_x_log10() +
      scale_y_log10() +
      theme_classic() +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
      theme(plot.title = element_text(hjust=0.5, face="bold"))
    return (data)

  }



}

scQCFilter = function(object,
                      minUMI = 500, minGenes = 250,
                      maxUMI = 100000, maxGenes = 50000,
                      minGenesPerUMI = 0.8 ,
                      maxGenesPerUMI = 1,
                      minMtRatio = 0,
                      maxMtRatio = 0.2,
                      project = 'Slice_culture',
                      removeZeros = FALSE){
  if (!require("Seurat", quietly = TRUE))
    install.packages("Seurat", library(Seurat) )

  #Get mt ratio
  object$percent.MT = mitoRatio(object)
  #Get novelty score
  object$log10genesPerUMI = nScore(object)

  object <- subset(x = object,
                   subset= (nCount_RNA >= minUMI) &
                     (nCount_RNA <= maxUMI) &
                     (nFeature_RNA >= minGenes) &
                     (nFeature_RNA <= maxGenes) &
                     (log10genesPerUMI > minGenesPerUMI) &
                     (log10genesPerUMI < maxGenesPerUMI) &
                     (percent.MT > minMtRatio) &
                     (percent.MT < maxMtRatio)
  )

  # Remove genes with zero counts
  # Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
  # Only keep genes expressed in more than 10 cells
  if (removeZeros == TRUE){
    counts <- GetAssayData(object = object, slot = "counts")
    nonzero <- counts > 0
    qcGene <- Matrix::rowSums(nonzero) >= 10
    counts <- counts[qcGene, ]

    object <- CreateSeuratObject(counts, meta.data = object@meta.data,
                                 project = project)
  }

  return (object)

}













