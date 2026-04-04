#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(ape)
library(ggsci)
library(ComplexHeatmap)
library(phytools)
library(phangorn)
library(grid)
library(gridExtra)
library(pheatmap)

mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)
cols <- c("steelblue", mycolors_tmp, "black")
#names(cols) <- c("Chimpanzee","AFR","AMR","EAS","EUR","SAS","internal")
names(cols) <- c("b","AFR","AMR","EAS","EUR","SAS","a")
#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")


#args = c("output/chr17/chr17-45585160-INV-706887/data/chr17-45585160-INV-706887.filter.noMissing.vcf.hap.wCMP.hap.tped", "output/chr17/chr17-45585160-INV-706887/data/chr17-45585160-INV-706887.filter.noMissing.vcf.hap.wCMP.hap.tfam", "chr17-45585160-INV-706887")


fin_tped <- args[1]
fin_hapID <- args[2]
INV_ID <- args[3]
outf_prefix <- args[4]

orig1_tped <- read.delim(fin_tped, header=F, sep='', na.strings = ".")
orig_tped <-  t(as.matrix(orig1_tped[,-c(1,2,3,4)]))
orig_tfam <- read.delim(fin_hapID, header=F, sep='')

orig_tfam$V7 <- apply(orig_tfam, 1, function(x) paste(unlist(strsplit(as.character(x),"_"))[2:3],collapse="_"))



## create a mapping between hapID and INV status
tmp_hap1 <- rep(paste(orig_tfam$V7, "1", sep='_'), each=2)
tmp_hap2 <- rep(paste(orig_tfam$V7, "2", sep='_'), each=2)

tmp_hap3 <- rep(paste(orig_tfam$V1, "1", sep='_'), each=2)
tmp_hap4 <- rep(paste(orig_tfam$V1, "2", sep='_'), each=2)

num_sample <- dim(orig_tfam)[1]
hapID <- rep("", num_sample)
hapID[ seq(1, num_sample*2, 2)] <- tmp_hap1[ seq(1, num_sample*2, 2)]
hapID[ seq(2, num_sample*2, 2)] <- tmp_hap2[ seq(2, num_sample*2, 2)]
hapID <- hapID[-2]

orig_hapID <- rep("", num_sample)
orig_hapID[ seq(1, num_sample*2, 2)] <- tmp_hap3[ seq(1, num_sample*2, 2)]
orig_hapID[ seq(2, num_sample*2, 2)] <- tmp_hap4[ seq(2, num_sample*2, 2)]
orig_hapID <- orig_hapID[-2]

index_INVsite <- which(grepl(INV_ID,orig1_tped$V2))
index_INVhap <- which(orig_tped[,index_INVsite] == 1)

mapping_hap_INV <- data.frame(hapID= hapID, INV=FALSE, orig_hapID= orig_hapID)
mapping_hap_INV$INV[index_INVhap] <- TRUE

row.names(orig_tped) <- orig_hapID


annoRow = data.frame(INV=ifelse(mapping_hap_INV[,2]==T, 1, 0))
row.names(annoRow) <- mapping_hap_INV$orig_hapID
p2 <- pheatmap(orig_tped, cluster_cols = F, clustering_distance_rows = "binary", clustering_method = "complete", color=c("darkblue","darkorange"), labels_col = "", labels_row = NULL, annotation_row =annoRow)


#outPDF = sprintf("unpolarized_treeHaplotype.%s.pdf", outf_prefix )
outPDF = sprintf("%s.polarizedHap.wCMP.pdf", outf_prefix )
ggsave(outPDF, p2, width=12, height = 9, useDingbats = F)
p2
dev.off()



