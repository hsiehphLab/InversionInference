#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggtree)
library(treeio)
library(ape)
library(ggsci)

mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)
cols <- c("gray", mycolors_tmp, "black")
names(cols) <- c("Chimpanzee","AFR","AMR","EAS","EUR","SAS","internal")
#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")

#args = c("output/chrX/chrX-72996081-INV-91407/data/chrX-72996081-INV-91407.filter.noMissing.vcf.hap.tped", "output/chrX/chrX-72996081-INV-91407/data/chrX-72996081-INV-91407.filter.noMissing.vcf.hap.tfam", "output/chrX/chrX-72996081-INV-91407/data/IDs_wPopID.hapID",  "output/chrX/chrX-72996081-INV-91407/iqtree/chrX-72996081-INV-91407.treefile",  "output/chrX/chrX-72996081-INV-91407/iqtree/chrX-72996081-INV-91407.timetree.nex", "chrX-72996081-INV-91407")

#args = c("output/chr8/chr8-2340295-INV-42058/data/chr8-2340295-INV-42058.filter.noMissing.vcf.hap.tped", "output/chr8/chr8-2340295-INV-42058/data/chr8-2340295-INV-42058.filter.noMissing.vcf.hap.tfam", "output/chr8/chr8-2340295-INV-42058/data/IDs_wPopID.hapID","output/chr8/chr8-2340295-INV-42058/iqtree/chr8-2340295-INV-42058.treefile", "output/chr8/chr8-2340295-INV-42058/iqtree/chr8-2340295-INV-42058.timetree.nex", "chr8-2340295-INV-42058")

#args = c("output/chr11/chr11-50136371-INV-206505/data/chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tped", "output/chr11/chr11-50136371-INV-206505/data/chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tfam", "output/chr11/chr11-50136371-INV-206505/data/IDs_wPopID.hapID", "output/chr11/chr11-50136371-INV-206505/iqtree/chr11-50136371-INV-206505.treefile", "output/chr1/chr1-13084312-INV-62181/iqtree/chr1-13084312-INV-62181.timetree.nex", "tmp")

#args = c("output/chr1/chr1-13084312-INV-62181/data/chr1-13084312-INV-62181.filter.noMissing.vcf.hap.tped", "output/chr1/chr1-13084312-INV-62181/data/chr1-13084312-INV-62181.filter.noMissing.vcf.hap.tfam", "output/chr1/chr1-13084312-INV-62181/data/IDs_wPopID.hapID", "output/chr1/chr1-13084312-INV-62181/iqtree/chr1-13084312-INV-62181.treefile", "output/chr1/chr1-13084312-INV-62181/iqtree/chr1-13084312-INV-62181.timetree.nex", "tmp")

#args = c("chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tped", "chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tfam", "IDs_wPopID.hapID", "../iqtree/chr11-50136371-INV-206505.treefile", "tmp")

fin_tped <- args[1]
fin_tfam <- args[2]
fin_hapID <- args[3]
fin_iqtree <- args[4]
fin_timetree <- args[5]
INVID <- args[6]
outf_prefix <- args[7]

#orig1_tped <- read.delim("../../PCA_hap/chr11_50136371_INV_206505.filter.noMissing.vcf.hap.tped", header=F, sep='')
#orig_tfam <- read.delim("../../PCA_hap/chr11_50136371_INV_206505.filter.noMissing.vcf.hap.tfam", header=F, sep='')

orig1_tped <- read.delim(fin_tped, header=F, sep='')
orig_tped <-  t(orig1_tped[,-c(1,2,3,4)])
orig_tfam <- read.delim(fin_tfam, header=F, sep='')

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

orig_hapID <- rep("", num_sample)
orig_hapID[ seq(1, num_sample*2, 2)] <- tmp_hap3[ seq(1, num_sample*2, 2)]
orig_hapID[ seq(2, num_sample*2, 2)] <- tmp_hap4[ seq(2, num_sample*2, 2)]

index_INVsite <- which(grepl(INVID,orig1_tped$V2))
index_INVhap <- which(orig_tped[,index_INVsite] == 1)

mapping_hap_INV <- data.frame(hapID= hapID, INV=FALSE, orig_hapID= orig_hapID)
mapping_hap_INV$INV[index_INVhap] <- TRUE
write.table(mapping_hap_INV, "../data/mapping_hap_INV.txt", quote = F, row.names = F, sep='\t')


#orig_hapID <- read.delim("IDs_wPopID.hapID", header=F, sep='')
orig_hapID <- read.delim(fin_hapID, header=F, sep='')
row.names(orig_tped) <- orig_hapID[,2]


#iqtree <- read.iqtree("chr11_50136371_INV_206505.treefile" )
# plot unrooted phylogeny
iqtree <- read.iqtree(fin_iqtree)

num_seqs <- length(iqtree@phylo$tip.label)

df_iqtree <- as.data.frame(as_tibble(iqtree))

Pop <- sapply(strsplit(df_iqtree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_iqtree)[1]-num_seqs))

df_iqtree$Pop <- Pop
df_iqtree$Pop <- ifelse(df_iqtree$Pop == "CMP", "Chimpanzee", df_iqtree$Pop)

df_iqtree$branchGroup <- ifelse(df_iqtree$label == "CMP", 2, 1)
df_iqtree$branchGroup <- as.factor(df_iqtree$branchGroup)
df_iqtree$INV <- mapping_hap_INV[ match(df_iqtree$label,mapping_hap_INV$orig_hapID), "INV"]
df_iqtree[df_iqtree$label == "CMP", "INV"] = FALSE

iqtree2 <- iqtree
iqtree@phylo$edge.length.mod <- iqtree@phylo$edge.length
iqtree@phylo$edge.length[1] = iqtree@phylo$edge.length.mod[1]*0.1

func_fort <- function(tr) tr %>% fortify
fort_tree <- func_fort(iqtree2@phylo)
upXlim <- max(fort_tree$x)
lowXlim <- min(fort_tree$x)

ggtree(iqtree, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_iqtree + geom_label2(aes(label=UFboot, subset=UFboot > 70), show.legend = FALSE) +
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + geom_treescale(offset = -2) +  
  scale_x_continuous(limits = c(lowXlim,upXlim/7.5)) +  
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + geom_tiplab(aes(label=label), linetype=NA, size=4, hjust=-0.2) +  guides(linetype="none", size="none", color=guide_legend(title="Population"), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm') )

ggsave(paste(outf_prefix,".treefile.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)



# plot rooted/dated iqtree
mytimetree <- read.beast(fin_timetree)

num_seqs <- length(mytimetree@phylo$tip.label)

df_mytimetree <- as.data.frame(as_tibble(mytimetree))

Pop <- sapply(strsplit(df_mytimetree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree)[1]-num_seqs))

df_mytimetree$Pop <- Pop
df_mytimetree$Pop <- ifelse(df_mytimetree$Pop == "CMP", "Chimpanzee", df_mytimetree$Pop)

df_mytimetree$branchGroup <- ifelse(df_mytimetree$Pop == "CMP", 2, 1)
df_mytimetree$branchGroup <- as.factor(df_mytimetree$branchGroup)
df_mytimetree$INV <- mapping_hap_INV[ match(df_mytimetree$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree[df_mytimetree$Pop == "CMP", "INV"] = FALSE
df_mytimetree$INV <- as.factor(df_mytimetree$INV)

p1 <- ggtree(mytimetree, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree + 
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17),breaks = levels(df_mytimetree$INV)) + geom_treescale(x=0.2) +
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_range("CI_date", color="gray",size=2,alpha=0.7) + 
#  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years") 

p2 <- revts(p1) 
p3 <- p2 + scale_x_continuous(breaks=c(-6:0), labels=abs(-6:0))



mytimetree2 <- treeio::drop.tip(mytimetree, "CMP")
num_seqs <- length(mytimetree2@phylo$tip.label)

df_mytimetree2 <- as.data.frame(as_tibble(mytimetree2))

Pop <- sapply(strsplit(df_mytimetree2[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree2)[1]-num_seqs))

df_mytimetree2$Pop <- Pop
df_mytimetree2$Pop <- ifelse(df_mytimetree2$Pop == "CMP", "Chimpanzee", df_mytimetree2$Pop)

df_mytimetree2$branchGroup <- ifelse(df_mytimetree2$Pop == "CMP", 2, 1)
df_mytimetree2$branchGroup <- as.factor(df_mytimetree2$branchGroup)
df_mytimetree2$INV <- mapping_hap_INV[ match(df_mytimetree2$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree2[df_mytimetree2$Pop == "CMP", "INV"] = FALSE
df_mytimetree2$INV <- as.factor(df_mytimetree2$INV) 


p4 <- ggtree(mytimetree2, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree2 + 
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17),breaks = levels(df_mytimetree2$INV)) + geom_treescale(x=0.1) + 
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_range("CI_date", color="gray",size=2,alpha=0.7) + 
  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree2$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV==T), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years")

upperTMRCA <- round(abs(min(unlist(mytimetree2@data$CI_date),na.rm = T)), digits = 1)
p5 <- revts(p4) 
breaks=seq(-upperTMRCA,0, upperTMRCA/2)
labels = abs(seq(-upperTMRCA,0,upperTMRCA/2))
for (i in 1:length(labels)) { if (i%%2==0) labels[i]=""}
p6 <- p5 + scale_x_continuous(breaks=breaks, labels=labels)

ggarrange(p3,p6,ncol=2,common.legend = T, widths = c(6,12))

ggsave(paste(outf_prefix,".timetree.iqtree.pdf",sep=''), height = 9, width = 18, limitsize=F, useDingbats=FALSE)
print (unlist(mytimetree2@data$CI_date))



# date tree using APE's penalized algorithm
rt_iqtree2 <- root(iqtree2@phylo, outgroup = "CMP")
numZeroBL = table(rt_iqtree2$edge.length == 0)[2]
if (numZeroBL > 0.3 * Ntip(rt_iqtree2)) {
#  rt_iqtree2$edge.length = ifelse(rt_iqtree2$edge.length == 0, min(min(rt_iqtree2$edge.length[rt_iqtree2$edge.length!=0])/1e4, 1e-10), rt_iqtree2$edge.length)
  rt_iqtree2$edge.length = rt_iqtree2$edge.length + 1e-10
}
my_node <- getMRCA(rt_iqtree2, tip=c("CMP", rt_iqtree2$tip.label[2])) 
mycalibration <- makeChronosCalib(phy = rt_iqtree2, node = my_node, age.min=4, age.max = 6 )
count_eval = 0
stopTrying = FALSE
list_mytimetree <- vector(mode="list", length=10)
for (i in 1:10){
  convergence = FALSE
  while (convergence != TRUE){
    count_eval = count_eval + 1
    print(count_eval)
    tryCatch(
      mytimetree <- chronos(rt_iqtree2, lambda = 1, calibration = mycalibration, model="discrete", control = chronos.control(epsilon = 1e-8, eval.max = 1e5, dual.iter.max = 1000, tol=1e-8, nb.rate.cat=1))
      ,error=function(e) e)
	if(inherits(mytimetree,"error")) next
    convergence = attr(mytimetree, "convergence")
    print(convergence)
    if (convergence == TRUE){
      list_mytimetree[[i]] <- mytimetree
      break
    }
    if (count_eval > 100){
      if (i!= 1){
        stopTrying = TRUE
        break
      }
    }
  }
  if (stopTrying == TRUE) break
}

list_mytimetree <- list_mytimetree[lapply(list_mytimetree,length)>0]
list_logLik <- unlist(lapply(list_mytimetree, function(x) attr(x,"PHII")$logLik))
mytimetree <- list_mytimetree[[which.max(list_logLik)]]

## plot dated (APE) tree 
write.tree(mytimetree, paste(outf_prefix,".mytimetree.nwk",sep=''))
mytimetree <- read.iqtree(paste(outf_prefix,".mytimetree.nwk",sep=''))


num_seqs <- length(mytimetree@phylo$tip.label)

df_mytimetree <- as.data.frame(as_tibble(mytimetree))

Pop <- sapply(strsplit(df_mytimetree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree)[1]-num_seqs))

df_mytimetree$Pop <- Pop
df_mytimetree$Pop <- ifelse(df_mytimetree$Pop == "CMP", "Chimpanzee", df_mytimetree$Pop)

df_mytimetree$branchGroup <- ifelse(df_mytimetree$label == "CMP", 2, 1)
df_mytimetree$branchGroup <- as.factor(df_mytimetree$branchGroup)
df_mytimetree$INV <- mapping_hap_INV[ match(df_mytimetree$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree[df_mytimetree$label == "CMP", "INV"] = FALSE


p1 <- ggtree(mytimetree, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree + 
  geom_label2(aes(label=UFboot, subset=UFboot > 70), show.legend = FALSE) +
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + geom_treescale(x=2) +
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years") 

p2 <- revts(p1) 
p3 <- p2 + scale_x_continuous(breaks=c(-6:0), labels=abs(-6:0))


#mytimetree2 <- drop.tip(mytimetree@phylo, "CMP")
keeptips <- mytimetree@phylo$tip.label[ mytimetree@phylo$tip.label != "CMP"]
mytimetree2 <- tree_subset(mytimetree, keeptips, levels_back = 100, group_node = F, root_edge = T)

num_seqs <- length(mytimetree2@phylo$tip.label)

df_mytimetree2 <- as.data.frame(as_tibble(mytimetree2))

Pop <- sapply(strsplit(df_mytimetree2[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_mytimetree2)[1]-num_seqs))

df_mytimetree2$Pop <- Pop
df_mytimetree2$Pop <- ifelse(df_mytimetree2$Pop == "CMP", "Chimpanzee", df_mytimetree2$Pop)

df_mytimetree2$branchGroup <- ifelse(df_mytimetree2$Pop == "CMP", 2, 1)
df_mytimetree2$branchGroup <- as.factor(df_mytimetree2$branchGroup)
df_mytimetree2$INV <- mapping_hap_INV[ match(df_mytimetree2$label,mapping_hap_INV$orig_hapID), "INV"]
df_mytimetree2[df_mytimetree2$Pop == "CMP", "INV"] = FALSE




p4 <- ggtree(mytimetree2, aes(color=Pop, linetype=branchGroup), size=1.2)  %<+% df_mytimetree2 + 
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + geom_treescale(x=.6) + 
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS")) + 
  geom_tiplab(aes(label=mapping_hap_INV[ match(df_mytimetree2$label, mapping_hap_INV$orig_hapID),"hapID"], subset=INV), linetype=NA, size=5, hjust=-0.1) +  
  guides(linetype="none", size="none", color=guide_legend(title="Population", nrow=1), shape=guide_legend(title=paste("Inversion",outf_prefix,sep="\n"), nrow=1)) + theme_tree2() + 
  theme(legend.position="top", legend.justification='left', legend.text = element_text(size=14), legend.title = element_text(size=12), legend.key.size = unit(1, 'cm'), axis.text.x = element_text(size=14)) + xlab("Million years")


p5 <- revts(p4) 
breaks=seq(-6,0,0.25)
labels = abs(seq(-6,0,0.25))
for (i in 1:length(labels)) { if (i%%2==0) labels[i]=""}
p6 <- p5 + scale_x_continuous(breaks=breaks, labels=labels)

ggarrange(p3,p6,ncol=2,common.legend = T, widths = c(12,6))

ggsave(paste(outf_prefix,".timetree.ape.pdf",sep=''), height = 9, width = 18, limitsize=F, useDingbats=FALSE)

