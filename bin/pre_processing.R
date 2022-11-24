#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)

# Create a parser
parser <- ArgumentParser("Parameters for Seurat object: Pre-processing")

parser$add_argument("--fileUmi", help = "Path to file UMI count matrix")
parser$add_argument("--fileHto", help = "Path to file HTO matrix")
parser$add_argument("--ndelim", help = "For the initial identity calss for each cell, delimiter for the cell's column name",
                    default = "_")
parser$add_argument("--rdsObject", help = "True if inputs are rds objects", default = FALSE)

parser$add_argument("--selectMethod", help = "Selection method", default = "mean.var.plot")
parser$add_argument("--numberFeatures", help = "Number of features to be used when finding variable features", 
                    type = integer, default = 2000)
parser$add_argument("--assay", help = "Choose assay between RNA or HTO",default = "HTO")

parser$add_argument("--normalisationMethod", help = "Normalisation method", default = "CLR")
parser$add_argument("--margin", help = "Margin for normalisation", type="numeric", default = 2)

parser$add_argument( "--OutputFile",help="Prefix of output files containing the output of HTODemux hashtag", default = "preprocessed")

args <- parser$parse_args()

#umi stands for the RNA matrix
if( args$rdsObject == TRUE){
  umi <- readRDS(args$fileUmi)
  counts <- readRDS(args$fileHto)
}else{
  umi <- Read10X(data.dir = args$fileUmi)
  counts <- Read10X(data.dir = args$fileHto)
}

#Identify which UMI corresponds to which hashtag.
joint.bcs <- intersect(colnames(umi), colnames(counts))
print(joint.bcs)

umi <- umi[, joint.bcs]
counts <- counts[, joint.bcs]

# Setup Seurat object
hashtag <- CreateSeuratObject(counts = umi, names.delim = args$ndelim)

# Add HTO data as a new assay independent from RNA
hashtag[[args$assay]] <- CreateAssayObject(counts = counts)

# Normalize HTO data
hashtag <- NormalizeData(hashtag, assay = args$assay, normalization.method = args$normalisationMethod, margin = args$margin)

#Save Results
saveRDS(hashtag, file = paste0(args$OutputFile, ".rds"))





