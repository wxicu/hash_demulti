#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for HTODemux Visualisation")
parser$add_argument("--HashtagPath",help="S4 object saved from the first part of HTODemux", default = NULL)
parser$add_argument("--assay",help="Name of the Hashtag assay HTO by default", default = "HTO")

#Output graphs - Ridge Plot
parser$add_argument("--ridgePlot", help = "Generates a ridge plot from the results, True to generate", default = FALSE)
parser$add_argument("--ridgeNCol", help = "Number of columns for ridgePlot", default = 3)

#Output graphs - Scatter Feature
parser$add_argument("--featureScatter",help = "Generates a ridge plot from the results, True to generate", default = FALSE)
parser$add_argument("--scatterFeat1", help = "Feature 1 for Feature Scatter Plot", default = NULL)
parser$add_argument("--scatterFeat2", help = "Feature 2 for Feature Scatter Plot", default = NULL)

#Output graphs - Violin Plot
parser$add_argument("--vlnPlot", help = "Generates a violin plot from the results, True to generate", default = FALSE)
parser$add_argument("--vlnFeatures", help = "Features to plot (gene expression, metrics, PC scores, anything that can be retreived by FetchData)", default = NULL)
parser$add_argument("--vlnLog", help = "plot the feature axis on log scale", default = TRUE)

#Output graphs - tSNE
parser$add_argument("--tSNE", help = "Generate a two dimensional tSNE embedding for HTOs", default = FALSE)
parser$add_argument("--tSNEIdents", help = "What should we remove from the object (we have Singlet,Doublet and Negative)", default = "Negative")
parser$add_argument("--tSNEInvert", help = "True or False", default = TRUE)
parser$add_argument("--tSNEVerbose", help = "True or False", default = FALSE)
parser$add_argument("--tSNEApprox", help = "True or False", default = FALSE)
parser$add_argument("--tSNEDimMax", help = "max -> number of donors ",type = "integer", default = 1)
parser$add_argument("--tSNEPerplexity", help = "value for perplexity", type = "integer",  default = 100)

#Output graphs - Heatmap
parser$add_argument("--heatMap", help = "Generate a Heatmap", default = FALSE)
parser$add_argument("--heatMapNcells", help ="value for number of cells", type = "integer",  default = 500)

args <- parser$parse_args()

hashtag <-readRDS(args$HashtagPath)

# Ridge Plot
# Group cells based on the max HTO signal
if(args$ridgePlot == TRUE){
  RidgePlot(hashtag, assay = args$assay, features = rownames(hashtag[[args$assay]]), ncol = args$ridgeNCol)
  ggsave("ridge.jpeg", device = 'jpeg', dpi = 500, height = 10, width = 10)
}

if(args$featureScatter == TRUE){
  FeatureScatter(hashtag, feature1 = args$scatterFeat1 , feature2 = args$scatterFeat2)
  ggsave('featureScatter.jpeg', device = 'jpeg',dpi = 500)
}

if(args$vlnPlot == TRUE){
  Idents(hashtag) <- "HTO_classification.global"
  VlnPlot(hashtag, features = args$vlnFeatures, pt.size = 0.1, log = args$vlnLog)
  ggsave('violinPlot.jpeg', device = 'jpeg',dpi = 500)
}

if(args$tSNE == TRUE){
  hashtag.subset <- subset(hashtag, idents = args$tSNEIdents, invert = args$tSNEInvert)
  DefaultAssay(hashtag.subset) <- args$assay
  hashtag.subset <- ScaleData(hashtag.subset, features = rownames(hashtag.subset),
                                   verbose = args$tSNEVerbose)
  hashtag.subset <- RunPCA(hashtag.subset, features = rownames(hashtag.subset), approx = args$tSNEApprox)
  hashtag.subset <- RunTSNE(hashtag.subset, dims = 1:args$tSNEDimMax, perplexity = args$tSNEPerplexity)
  DimPlot(hashtag.subset)
  ggsave('tSNE.jpeg', device = 'jpeg',dpi = 500)
}

if(args$heatMap == TRUE){
  HTOHeatmap(hashtag, assay = args$assay, ncells = args$heatMapNcells)
  ggsave('heatMap.jpeg', device = 'jpeg',dpi = 500)
}





