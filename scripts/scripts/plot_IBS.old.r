#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(ggpubr)
library(ggsci)
library(reshape2)
library(ComplexHeatmap)
#library(dplyr)
#library(tidyverse)
#library(rstatix)   
#source("~/bin/scripts/scripts_uw/R_plotting/R_rainclouds.r")
source("/Users/hseihph/mnt/vol26/home/hsiehph/bin/scripts/scripts_uw/R_plotting/R_rainclouds.r")

#ibs <- read.delim("tmp.out")
#elemDiff_pairHaps <- read.delim("chr11-50136371-INV-206505.elemDiff_pairHaps.dat", header=F, sep='')
ibs <- read.delim(args[1])
elemDiff_pairHaps <- read.delim(args[2], header=F , sep='')
outPDF1 <- args[3]
outPDF2 <- args[4]

ibs$INVgeno <- paste(ibs$h1_INVgeno, ibs$h2_INVgeno, sep="")
ibs$INVgeno <- ifelse(ibs$INVgeno == "10", "01", ibs$INVgeno)

min_ibs = round(min(ibs$IBS), digits=1) - 0.05
max_ibs = round(max(ibs$IBS), digits=1) + 0.05
tick_breaks = seq(max(0,round(min_ibs,1)), min(round(max_ibs,1),1), 0.1)

p1 <- ggplot(ibs, aes(x=INVgeno, y=IBS, fill=INVgeno)) + coord_flip() + 
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .6, adjust =1, trim = T) +
  geom_point(position=position_jitter(width = .15), size = 1, alpha = 0.5, shape=21) +
  geom_boxplot(aes(x=INVgeno, y=IBS), outlier.shape = NA, alpha=0.7, width=0.1, color="black") +
  scale_color_nejm() + scale_fill_nejm() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank(), panel.border = element_blank(), 
        axis.line = element_line(colour = "black"), 
        axis.title=element_text(size=16), axis.text=element_text(size=16)) +
  ylab("Pairwise identity by state") + xlab("Inversion genotype") + guides(fill=FALSE) + 
  scale_y_continuous(breaks=tick_breaks, limits = c(min_ibs,max_ibs))


pdf(outPDF1, height = 6, width = 9)
#plot_grid(p1, p2, ncol=1, align = "v", rel_heights = c(3,1))
p1
dev.off()


# heatmap element-wise difference between pairs of haplotypes
cols = pal_jama()(2)
colnames = c("False","True")


elemDiff_pairHaps <- elemDiff_pairHaps[ order(elemDiff_pairHaps$V1, elemDiff_pairHaps$V2),]

row.names(elemDiff_pairHaps) <- paste(elemDiff_pairHaps$V3, elemDiff_pairHaps$V4,sep=",")

list_inv_anno <- paste(elemDiff_pairHaps$V1, elemDiff_pairHaps$V2, sep="")
list_inv_anno <- ifelse(list_inv_anno=="10", "01", list_inv_anno)

start_groups <- c(1, unname(cumsum(table(list_inv_anno)))+1)[-4]
end_groups <- unname(cumsum(table(list_inv_anno)))

m3 = c("00","01","11")
m3_labels = c("D/D","D/I","I/I")

m_elemDiff_pairHaps <- as.matrix(elemDiff_pairHaps[,-(1:4)])

ht_list <- NULL
for (i in 1:3){
  row_split = c(rep(1,sum(list_inv_anno==m3[i])))
  tmp_m_elemDiff_pairHaps <- m_elemDiff_pairHaps[start_groups[i]:end_groups[i],]
  anno_right = rowAnnotation(foo=anno_block(gp=gpar(fill=pal_d3()(3)[i]), labels = m3_labels[i], labels_gp = gpar(col="white", fontsize=10)))

  if(i==2){rowTitle = "Haplotype pairs"} else{rowTitle=""}
  ht_list <- c(ht_list, Heatmap(tmp_m_elemDiff_pairHaps, cluster_rows = T, cluster_columns = F, name="Mismatch", col =  cols, heatmap_legend_param = list(at = c("0", "1"), labels = colnames), show_column_names = F, show_row_names = F, column_title = "SNV", row_title = rowTitle, right_annotation = anno_right))
}

all_ht_list <- ht_list[[1]] %v% ht_list[[2]] %v% ht_list[[3]]
pdf(outPDF2, height = 6, width = 9)
all_ht_list
dev.off()
