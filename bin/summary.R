#!/usr/bin/env Rscript
library(data.table)
library(stringr)
library(argparse)
library(dplyr)

parser <- ArgumentParser("Parameters for comparing parameters")
parser$add_argument("--demuxem", help = "Folder containing output files of demuxem", default = NULL)
parser$add_argument("--htodemux", help = "Folder containing output files of htodemux", default = NULL)
parser$add_argument("--multiseq", help = "Folder containing output files of multiseq", default = NULL)
parser$add_argument("--hashsolo", help = "Folder containing output files of hashsolo", default = NULL)
parser$add_argument("--solo", help = "Folder containing output files of solo", default = NULL)
parser$add_argument("--hashedDrops", help = "Folder containing output files of hashedDrops", default = NULL)

args <- parser$parse_args()

demuxem_summary <- function(demuxem_res) {
  assign <- lapply(demuxem_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_obs.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    colnames(obs_res)[1] <- "Barcode"
    demuxem_assign <- obs_res[, c("Barcode", "assignment")]
    colnames(demuxem_assign)[2] <- basename(x)
    demuxem_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "demuxem_assignment.csv", row.names=FALSE)
  
  classi <- lapply(demuxem_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_obs.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    colnames(obs_res)[1] <- "Barcode"
    demuxem_classi <- obs_res[, c("Barcode", "demux_type")]
    colnames(demuxem_classi)[2] <- basename(x)
    demuxem_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "demuxem_classification.csv", row.names=FALSE)
  
  params <- lapply(demuxem_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "demuxem_params.csv", row.names=FALSE)
  
}

hashsolo_summary <- function(hashsolo_res){
  assign <- lapply(hashsolo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    hashsolo_assign <- obs_res[, c("V1", "Classification")]
    colnames(hashsolo_assign) <- c("Barcode", basename(x))
    hashsolo_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "hashsolo_assignment.csv", row.names=FALSE)
  
  classi <- lapply(hashsolo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    hashsolo_classi <- obs_res[, c("V1", "most_likely_hypothesis")]
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 0]  <- "negative"
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 1]  <- "singlet"
    hashsolo_classi$most_likely_hypothesis[hashsolo_classi$most_likely_hypothesis == 2]  <- "doublet"
    colnames(hashsolo_classi) <- c("Barcode", basename(x))
    hashsolo_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "hashsolo_classification.csv", row.names=FALSE)
  
  params <- lapply(hashsolo_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "hashsolo_params.csv", row.names=FALSE)
}

hasheddrops_summary <- function(hasheddrops_res){
  classi <- lapply(hasheddrops_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    obs_res <- fread(obs_res_dir, header = TRUE)
    obs_res$Classification = ifelse(obs_res$Confident, "Singlet", ifelse(obs_res$Doublet, "Doublet", "Negative"))
    colnames(obs_res)[1] <- 'Barcode'
    hasheddrops_classi <- obs_res[,c('Barcode', 'Classification')]
    colnames(hasheddrops_classi) <- c("Barcode", basename(x))
    hasheddrops_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "hasheddrops_classification.csv", row.names=FALSE)
  
  params <- lapply(hasheddrops_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "hasheddrops_params.csv", row.names=FALSE)
}

multiseq_summary <- function(multiseq_res){
  assign <- lapply(multiseq_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    multiseq_assign <- fread(obs_res_dir, header = TRUE)
    colnames(multiseq_assign) = c("Barcode", basename(x))
    multiseq_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "multiseq_assignment.csv", row.names=FALSE)
  
  params <- lapply(multiseq_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "multiseq_params.csv", row.names=FALSE)
}

htodemux_summary <- function(htodemux_res){
  assign <- lapply(htodemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_assignment_htodemux.csv", full.names = TRUE)[1]
    htodemux_assign <- fread(obs_res_dir, header = TRUE)
    colnames(htodemux_assign) = c("Barcode", basename(x))
    htodemux_assign
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(assign, "htodemux_assignment.csv", row.names=FALSE)
  
  classi <- lapply(htodemux_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_classification_htodemux.csv", full.names = TRUE)[1]
    htodemux_classi <- fread(obs_res_dir, header = TRUE)
    colnames(htodemux_classi) = c("Barcode", basename(x))
    htodemux_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "htodemux_classification.csv", row.names=FALSE)
  
  params <- lapply(htodemux_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    params_res <- params_res[,2:3]
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "htodemux_params.csv", row.names=FALSE)
}

solo_summary <- function(solo_res){
  classi <- lapply(solo_res, function(x){
    obs_res_dir <- list.files(x, pattern = "_res.csv", full.names = TRUE)[1]
    solo_classi <- fread(obs_res_dir, header = TRUE)
    colnames(solo_classi) = c("Barcode", basename(x))
    solo_classi$Barcode <- gsub("-0", "", solo_classi$Barcode)
    solo_classi
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
  write.csv(classi, "solo_classification.csv", row.names=FALSE)
  
  params <- lapply(solo_res, function(x){
    params_dir <- list.files(x, pattern = "params.csv", full.names = TRUE)[1]
    params_res <- fread(params_dir, header = TRUE)
    colnames(params_res)[2] <- basename(x)
    params_res
  }) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Argument"), .)
  write.csv(params, "solo_params.csv", row.names=FALSE)
}


if (!is.null(args$hashedDrops)){
  hashedDrops_res <- substring(args$hashedDrops, 1, nchar(args$hashedDrops)-1)
  hashedDrops_res <- str_split(hashedDrops_res, pattern=':')[[1]]
  hasheddrops_summary(hashedDrops_res)
  print("hashedDrops result found")
}
if (!is.null(args$demuxem)){
  demuxem_res <- substring(args$demuxem, 1, nchar(args$demuxem)-1)
  demuxem_res <- str_split(demuxem_res, pattern=':')[[1]]
  demuxem_summary(demuxem_res)
  print("DemuxEM result found")
}
if (!is.null(args$hashsolo)){
  hashsolo_res <- substring(args$hashsolo, 1, nchar(args$hashsolo)-1)
  hashsolo_res <- str_split(hashsolo_res, pattern=':')[[1]]
  hashsolo_summary(hashsolo_res)
  print("HashSolo result found")
}
if (!is.null(args$multiseq)){
  multiseq_res <- substring(args$multiseq, 1, nchar(args$multiseq)-1)
  multiseq_res <- str_split(multiseq_res, pattern=':')[[1]]
  multiseq_summary(multiseq_res)
  print("MultiSeqDemux result found")
}
if (!is.null(args$htodemux)){
  htodemux_res <- substring(args$htodemux, 1, nchar(args$htodemux)-1)
  htodemux_res <- str_split(htodemux_res, pattern=':')[[1]]
  htodemux_summary(htodemux_res)
  print("HTODemux result found")
}
if (!is.null(args$solo)){
  solo_res <- substring(args$solo, 1, nchar(args$solo)-1)
  solo_res <- str_split(solo_res, pattern=':')[[1]]
  solo_summary(solo_res)
  print("solo result found")
}

assignment <- list.files(".", pattern = "_assignment.csv", full.names = TRUE)
assignment_all <- lapply(assignment, function(x){
  assign <- fread(x, header = TRUE)
  assign
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
write.csv(assignment_all, "assignment_all.csv", row.names=FALSE)

classification <- list.files(".", pattern = "_classification.csv", full.names = TRUE)
classification_all <- lapply(classification, function(x){
  classi <- fread(x, header = TRUE)
  classi
}) %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Barcode"), .)
write.csv(classification_all, "classification_all.csv", row.names=FALSE)

# demuxem_res <- str_split("../hash_out/demuxem/demuxem_1;../hash_out/demuxem/demuxem_2", pattern=';')[[1]]
# demuxem_summary(demuxem_res)
# hashsolo_res <- str_split("../hash_out/hashsolo/hashsolo_1;../hash_out/hashsolo/hashsolo_2", pattern=';')[[1]]
# hashsolo_summary(hashsolo_res)
# hasheddrops_res <- str_split("../hash_out/hashedDrops/hashedDrops_1;../hash_out/hashedDrops/hashedDrops_2", pattern=';')[[1]]
# hasheddrops_summary(hasheddrops_res)
# multiseq_res <- str_split("../hash_out/multiseq/multiseq_1;../hash_out/multiseq/multiseq_2", pattern=';')[[1]]
# multiseq_summary(multiseq_res)
# htodemux_res <- str_split("../hash_out/htodemux/htodemux_1;../hash_out/htodemux/htodemux_2", pattern=';')[[1]]
# htodemux_summary(htodemux_res)

# solo <- fread("~/Desktop/hash_out/solo/solo_predict.csv", header = TRUE)
# colnames(solo) <- c("Barcode", "Classification")
# solo$Barcode <- gsub("*-0","",solo$Barcode)
# 
# htodemul_assign = fread("../hash_out/htodemux/htodemux_out/htodemux_assignment_htodemux.csv", header = TRUE) 
# colnames(htodemul_assign) = c("Barcode","Assignment")
# htodemul_assign
# 
# htodemul_classi = fread("../hash_out/htodemux/htodemux_out/htodemux_classification_htodemux.csv", header = TRUE)  
# colnames(htodemul_classi) = c("Barcode", "Classification")
# htodemul_classi
# 
# multiseq_assign = fread("../hash_out/multiseq/multiseq_out/multiseq.csv", header = TRUE)  
# colnames(multiseq_assign) = c("Barcode","Assignment")
# multiseq_assign
# 
# hasheddrops_res = fread("../hash_out/hashedDrops/hashedDrops_out/hashedDrops.csv", header = TRUE)  
# hasheddrops_res$Classification = ifelse(hasheddrops_res$Confident, "Singlet", ifelse(hasheddrops_res$Doublet, "Doublet", "Negative"))
# colnames(hasheddrops_res)[1] = 'Barcode'
# hasheddrops_classi = hasheddrops_res[,c('Barcode', 'Classification')]
# hasheddrops_classi
# 
# demuxem_res = fread("../hash_out/demuxem/demuxem_out/demuxem_res_obs.csv", header = TRUE)  
# colnames(demuxem_res)[1] = "Barcode"
# demuxem_assign <- demuxem_res[, c("Barcode", "assignment")]
# demuxem_assign
# demuxem_classi <- demuxem_res[, c("Barcode", "demux_type")]
# colnames(demuxem_classi)[2] = "Classification"
# demuxem_classi
# 
# hash_solo <- fread("~/Desktop/hash_out/hashsolo/hashsolo_out/hashsolo.csv")
# hash_solo <- hash_solo[, c("V1", "most_likely_hypothesis", "Classification")]
# colnames(hash_solo) <- c("Barcode", "Classification", "Assignment")
# hash_solo$Classification[hash_solo$Classification == 0]  <- "negative"
# hash_solo$Classification[hash_solo$Classification == 1]  <- "singlet"
# hash_solo$Classification[hash_solo$Classification == 2]  <- "doublet"
# hash_solo_assign <- hash_solo[, -c("Classification")]
# hash_solo_classi <- hash_solo[, -c("Assignment")]
# table(hash_solo$Assignment)

# a="/Users/xichenwu/hash_demulti的副本/work/da/cf51271552727c547df0897b702418/multiseq_1;"
# a = substring(a,1, nchar(a)-1)
# multiseq_res <- str_split(a, pattern=';')[[1]]
# multiseq_summary(multiseq_res)

#a = "/Users/xichenwu/hash_demulti/hash_out/solo/solo_1"
#solo_res <- str_split(a, pattern=';')[[1]]
#solo_summary(solo_res)
