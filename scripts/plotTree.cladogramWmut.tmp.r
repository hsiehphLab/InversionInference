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

args = c("temp.nwk", "../candidates_20kbAND2Mbp_auto_v3/output/chr11/chr11-50136371-INV-206505/data/mapping_hap_INV.txt", "11p11.4.treefile.wMut.nwk.pdf.cladogramMut.pdf", "output/chr11/chr11-50136371-INV-206505/data/df_hclust.txt")


args = c("chr3-195749464-INV-230745.treefile.wMut.nwk", "../candidates_20kbAND2Mbp_auto_v3/output/chr3/chr3-195749464-INV-230745/data/mapping_hap_INV.txt", "chr3-195749464-INV-230745.treefile.wMut.nwk.pdf.cladogramMut.pdf", "output/chr3/chr3-195749464-INV-230745/data/df_hclust.txt")


#args = c("../candidates_above2Mb_auto_100kbChunks_chr8only/output/chr8/chr8-7301025-INV-5297356_9/iqtree_btFASTA/bt.55.timetree.nwk", "../candidates_above2Mb_auto_100kbChunks_chr8only/output/chr8/chr8-7301025-INV-5297356_9/data/mapping_hap_INV.txt","8p23.1.9.pdf.cladogramMut.pdf")


fin_tree <- args[1]
fin_mapping_hap_INV <- args[2]
outf_prefix <- args[3]

fin_df_hclust <- args[4]
df_hclust <- read.delim(fin_df_hclust, sep='', header = F)

mapping_hap_INV <- read.delim(fin_mapping_hap_INV, sep='')

tree <- treeio::drop.tip(read.tree(fin_tree), "CMP")

df_tree <- as.data.frame(as_tibble(tree))

getState <- function(r){
  a <- unlist(strsplit(r[4], split="_"))
  if (length(a) > 2){
    return(NA)
  }else{
    return(a[2])
  }
}
num_seqs = length(tree$tip.label)
df_tree$state <- apply(df_tree, 1, function(x) getState(x))
df_tree$numState <- apply(df_tree,1, function(x) ifelse(grepl("notMut",x[5]), NA, nchar(x[5])))
df_tree$invEvent <- ifelse(df_tree$numState>1, as.character(1), NA)

df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]
df_tree$INVcode <- ifelse(df_tree$INV == T, as.character(2), NA)
df_tree$INVcolor <- df_hclust[ match(df_tree$label, df_hclust$V3), "V4"]

mycolors <- c("#7F8080", "#2081f9",pal_jama()(7))
names(mycolors) <- c("branch","invEvent",1:7)  
  
ggtree(tree, layout = "circular", branch.length = "none", aes(color="branch"), size=1.2)  %<+% df_tree + 
  geom_tippoint(aes(shape=INV, color=as.factor(INVcolor)), size=5) +
  geom_text(aes(label=INVcolor)) + 
  scale_shape_manual(values=c(20,6,17)) + 
  geom_nodepoint(aes(color="invEvent", shape=invEvent), size=8) +
#  scale_color_manual(values=c("#7F8080","#2081f9", "orange","#A51E2B")) + 
  scale_color_manual(values=mycolors) + 
  guides(fill="none", size="none", color=guide_legend(title="Population"), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix,".cladogramMut.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)

