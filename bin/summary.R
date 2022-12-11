library(data.table)
library(stringr)
a = TRUE
list.files("hash_demulti的副本", pattern = "\\.nf$", full.names = TRUE)[1]
hash_solo <- fread("~/Desktop/hash_out/hashsolo/hashsolo_out/hashsolo.csv")
hash_solo <- hash_solo[, c("V1", "most_likely_hypothesis", "Classification")]
colnames(hash_solo) <- c("Barcode", "Classification", "Assignment")
hash_solo$Classification[hash_solo$Classification == 0]  <- "negative"
hash_solo$Classification[hash_solo$Classification == 1]  <- "singlet"
hash_solo$Classification[hash_solo$Classification == 2]  <- "doublet"
hash_solo_assign <- hash_solo[, -c("Classification")]
hash_solo_classi <- hash_solo[, -c("Assignment")]
table(hash_solo$Assignment)

solo <- fread("~/Desktop/hash_out/solo/solo_predict.csv", header = TRUE)
colnames(solo) <- c("Barcode", "Classification")
solo$Barcode <- gsub("*-0","",solo$Barcode)
res <- merge(solo, hash_solo, all=TRUE)

htodemul_assign = fread("../hash_out/htodemux/htodemux_out/htodemux_assignment_htodemux.csv", header = TRUE) 
colnames(htodemul_assign) = c("Barcode","Assignment")
htodemul_assign

htodemul_classi = fread("../hash_out/htodemux/htodemux_out/htodemux_classification_htodemux.csv", header = TRUE)  
colnames(htodemul_classi) = c("Barcode", "Classification")
htodemul_classi

multiseq_assign = fread("../hash_out/multiseq/multiseq_out/multiseq.csv", header = TRUE)  
colnames(multiseq_assign) = c("Barcode","Assignment")
multiseq_assign

hasheddrops_res = fread("../hash_out/hashedDrops/hashedDrops_out/hashedDrops.csv", header = TRUE)  
hasheddrops_res$Classification = ifelse(hasheddrops_res$Confident, "Singlet", ifelse(hasheddrops_res$Doublet, "Doublet", "Negative"))
colnames(hasheddrops_res)[1] = 'Barcode'
hasheddrops_classi = hasheddrops_res[,c('Barcode', 'Classification')]
hasheddrops_classi

demuxem_res = fread("../hash_out/demuxem/demuxem_out/demuxem_res_obs.csv", header = TRUE)  
colnames(demuxem_res)[1] = "Barcode"
demuxem_assign <- demuxem_res[, c("Barcode", "assignment")]
demuxem_assign
demuxem_classi <- demuxem_res[, c("Barcode", "demux_type")]
colnames(demuxem_classi)[2] = "Classification"
demuxem_classi
