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


args = c("8p23_left/8201024-8301024.nex.nwk", "8p23_left/mapping_hap_INV.txt", "chr8-7301025-INV-5297356_distal")

fin_tree <- args[1]
fin_mapping_hap_INV <- args[2]
outf_prefix <- args[3]

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
df_tree$numState <- apply(df_tree,1, function(x) nchar(x[5]))
df_tree$invEvent <- ifelse(df_tree$numState>1, as.character(1), NA)

df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]
df_tree$INVcode <- ifelse(df_tree$INV == T, as.character(2), NA)


ggtree(tree, layout = "circular", branch.length = "none", aes(color="#7F8080"), size=1.2)  %<+% df_tree + 
  geom_tippoint(aes(shape=INV, color=INV), size=5) + 
  scale_shape_manual(values=c(20,NA,17,NA)) + 
  geom_nodepoint(aes(color=invEvent, shape=invEvent), size=8) +
  scale_color_manual(values=c("#7F8080","#2081f9", "#7F8080","#A51E2B")) + 
  guides(fill="none", size="none", color=guide_legend(title="Population"), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix,".distalTree.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)


## proximal 

args = c("8p23_left/12001024-12101024.nex.nwk", "8p23_left/mapping_hap_INV.txt", "chr8-7301025-INV-5297356_proximal")

fin_tree <- args[1]
fin_mapping_hap_INV <- args[2]
outf_prefix <- args[3]

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
df_tree$numState <- apply(df_tree,1, function(x) nchar(x[5]))
df_tree$invEvent <- ifelse(df_tree$numState>1, as.character(1), NA)

df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]
df_tree$INVcode <- ifelse(df_tree$INV == T, as.character(2), NA)


ggtree(tree, layout = "circular", branch.length = "none", aes(color="#7F8080"), size=1.2)  %<+% df_tree + 
  geom_tippoint(aes(shape=INV, color=INV), size=5) + 
  scale_shape_manual(values=c(20,NA,17,NA)) + 
  geom_nodepoint(aes(color=invEvent, shape=invEvent), size=8) +
  scale_color_manual(values=c("#7F8080","#2081f9", "#7F8080","#A51E2B")) + 
  guides(fill="none", size="none", color=guide_legend(title="Population"), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix,".proximalTree.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)


args = c("8p23_left/chr8-7301025-INV-5297356.nex.nwk", "8p23_left/mapping_hap_INV.txt", "chr8-7301025-INV-5297356")

fin_tree <- args[1]
fin_mapping_hap_INV <- args[2]
outf_prefix <- args[3]

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
df_tree$numState <- apply(df_tree,1, function(x) nchar(x[5]))
df_tree$invEvent <- ifelse(df_tree$numState>1, as.character(1), NA)

df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]
df_tree$INVcode <- ifelse(df_tree$INV == T, as.character(2), NA)


ggtree(tree, layout = "circular", branch.length = "none", aes(color="#7F8080"), size=1.2)  %<+% df_tree + 
  geom_tippoint(aes(shape=INV, color=INV), size=5) + 
  scale_shape_manual(values=c(20,NA,17,NA)) + 
  geom_nodepoint(aes(color=invEvent, shape=invEvent), size=8) +
  scale_color_manual(values=c("#7F8080","#2081f9", "#7F8080","#A51E2B")) + 
  guides(fill="none", size="none", color=guide_legend(title="Population"), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix,".wholeTree.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)






## plot karyotypes.
library(karyoploteR)
plot.params <- getDefaultPlotParams(plot.type=1)
plot.params$ideogramheight <- 10

quartz(type="pdf", file="chr8p23_karyoty.pdf", height=9, width=12)
chromosomes <- "chr8"
kp <- plotKaryotype(plot.type = 1, chromosomes=chromosomes, plot.params = plot.params)
kpPlotRegions(kp, data=c("chr8:8201024-8301024", "chr8:12001024-12101024"), col="#AACCFF", r0=0, r1=0.25)
dev.off()

quartz(type="pdf", file="chr17q21_karyoty.pdf", height=9, width=12)
chromosomes <- "chr17"
kp <- plotKaryotype(plot.type = 1, chromosomes=chromosomes, plot.params = plot.params)
kpPlotRegions(kp, data=c("chr17:45585159-46292045"), col="#AACCFF", r0=0, r1=0.25)
dev.off()


quartz(type="pdf", file="chr11q11_karyoty.pdf", height=9, width=12)
chromosomes <- "chr11"
kp <- plotKaryotype(plot.type = 1, chromosomes=chromosomes, plot.params = plot.params)
kpPlotRegions(kp, data=c("chr11:50154999-50324102"), col="#AACCFF", r0=0, r1=0.25)
dev.off()



