#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for HTODemux")
parser$add_argument("--seuratObjectPath", help = "seurat object")
parser$add_argument("--quantile", help = "Positive quantile per default: 0.99", type = "double", default = 0.99)
parser$add_argument("--kfunc", help = "Cluster function choose between: Clara - kmeans", type = "character", default = "clara")
parser$add_argument("--nstarts", help = "number of starts for demultiplex", type = "integer", default = 100)
parser$add_argument("--nsamples", help = "number of samples for demultiplex", type = "integer", default = 100)
parser$add_argument("--seed", help = "sets random seed", type = "integer", default = 42)
parser$add_argument("--init", help = "Initial number of clusters for hashtags", default = NULL)
parser$add_argument("--objectOutHTO", help = "Name for the object containing the output of HTODemux object", default = "result")
parser$add_argument("--assignmentOutHTO", help="Name for the file containing the output of HTODemux assignment", default = "result")

args <- parser$parse_args()

# Loading Seurat object
hashtag <-readRDS(args$seuratObjectPath)

# Demultiplex cells based on HTO enrichment
hashtag <- HTODemux(hashtag, assay = "HTO", positive.quantile = args$quantile, nstarts = args$nstarts, kfunc = args$kfunc)

# Global classification results
table(hashtag$HTO_classification.global)


table(hashtag$HTO_classification)


print("-----------------------------------------------")
#dim(x = hashtag)
#head(x = rownames(x = hashtag))
#head(x = colnames(x = hashtag))
#names(x = hashtag)
#colnames(x = hashtag[[]])

hashtag[['RNA']]
print("-----------------------------------------------")
hashtag[['HTO']]

print("-----------------------------------------------")


#print("------------------- Percentage of largest gene ----------------------------")
hashtag

# Saving results

print("------------------- Following Files are saved ----------------------------")
print(paste0(args$assignmentOutHTO, "_assignment_htodemux.csv"))
print(paste0(args$assignmentOutHTO, "_classification_htodemux.csv"))
print(paste0(args$objectOutHTO,".rds"))
write.csv(hashtag$HTO_classification, paste0(args$assignmentOutHTO, "_assignment_htodemux.csv"))
write.csv(hashtag$HTO_classification.global, paste0(args$assignmentOutHTO, "_classification_htodemux.csv"))
saveRDS(hashtag, file=paste0(args$objectOutHTO,".rds"))


