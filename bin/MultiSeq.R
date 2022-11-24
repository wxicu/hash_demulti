#!/usr/bin/env Rscript
#Libraries
library(Seurat)
library(argparse)
library(ggplot2)

# Create a parser
parser <- ArgumentParser("Parameters for MultiSeqDemux")
parser$add_argument("--seuratObjectPath", help = "Seurat object")
parser$add_argument("--quantile", help = "The quantile to use for classification", type = "double", default = 0.7)
parser$add_argument("--autoThresh", help = "Whether to perform automated threshold finding to define the best quantile", default = TRUE)
parser$add_argument("--maxiter", help = "Maximum number of iterations if autoThresh = TRUE ", type = "double", default = 5)
parser$add_argument("--qrangeFrom", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.1)
parser$add_argument("--qrangeTo", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.9)
parser$add_argument("--qrangeBy", help = "A range of possible quantile values to try if autoThresh is TRUE", type = "double", default = 0.05)
parser$add_argument("--verbose", help = "Prints the output", default = TRUE)
parser$add_argument("--assay", help = "Name of the multiplexing assay (HTO by default)")
parser$add_argument("--assignmentOutMulti",help = "Name for the file containing the output of MULTI-Seq Demux object", default = "result.csv")
parser$add_argument("--objectOutMulti", help = "Name for the object containing the output of MULTI-Seq Demux object", default = "result")
args <- parser$parse_args()
hashtag <-readRDS(args$seuratObject)
if (args$autoThresh == TRUE) {
    hashtag <- MULTIseqDemux(hashtag, assay = args$assay,  quantile = args$quantile, autoThresh = TRUE, maxiter = args$maxiter, qrange=seq(from = args$qrangeFrom, to = args$qrangeTo, by = args$qrangeBy), verbose = args$verbose)
}else{
    hashtag <- MULTIseqDemux(hashtag, assay = args$assay,  quantile = args$quantile, qrange=seq(from = args$qrangeFrom, to = args$qrangeTo, by = args$qrangeBy), verbose = args$verbose)
    
}

table(hashtag$MULTI_ID)
print("----------------------------------------------------------------------------")
table(hashtag$MULTI_classification)

hashtag
print("-----------------------------------------------")

#dim(x = hashtag)
head(x = rownames(x = hashtag))

head(x = colnames(x = hashtag))

names(x = hashtag)

print("-----------------------------------------------")

hashtag[['RNA']]
print("-----------------------------------------------")
hashtag[['HTO']]

print("-----------------------------------------------")

colnames(x = hashtag[[]])



#Save Results
print("------------------- Following Files are saved ----------------------------")
print(paste0(args$assignmentOutMulti, "_assignment_multiseq_demux.csv"))
print(paste0(args$objectOutMulti",".rds"))
write.csv(hashtag$MULTI_ID, paste0(args$assignmentOutMulti,".csv"))
saveRDS(hashtag, paste0(args$objectOutMulti",".rds"))
