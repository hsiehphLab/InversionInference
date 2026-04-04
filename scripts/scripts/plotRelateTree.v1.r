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

mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)
cols <- c("steelblue", mycolors_tmp, "black")
#names(cols) <- c("Chimpanzee","AFR","AMR","EAS","EUR","SAS","internal")
names(cols) <- c("b","AFR","AMR","EAS","EUR","SAS","a")
#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")


#args = c("output/chr11/chr11-50136371-INV-206505/relate/chr11-50136371-INV-206505.filter.haps", "output/chr11/chr11-50136371-INV-206505/relate/chr11-50136371-INV-206505.filter.poplabels", "chr11-50136371-INV-206505", "output/chr11/chr11-50136371-INV-206505/relate/all.nwk", "tmp")

#args = c("output/chr17/chr17-45585160-INV-706887/relate/chr17-45585160-INV-706887.filter.haps", "output/chr17/chr17-45585160-INV-706887/relate/chr17-45585160-INV-706887.filter.poplabels", "chr17-45585160-INV-706887", "output/chr17/chr17-45585160-INV-706887/relate/all.nwk", "tmp")

#args = c("output/chr1/chr1-13084312-INV-62181/relate/chr1-13084312-INV-62181.filter.haps", "output/chr1/chr1-13084312-INV-62181/relate/chr1-13084312-INV-62181.filter.poplabels", "chr1-13084312-INV-62181", "output/chr1/chr1-13084312-INV-62181/relate/all.nwk", "tmp")

#args = c("chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tped", "chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tfam", "IDs_wPopID.hapID", "../iqtree/chr11-50136371-INV-206505.treefile", "tmp")

#args = c("output/chr8/chr8-7301025-INV-5297356/relate/chr8-7301025-INV-5297356.filter.haps","output/chr8/chr8-7301025-INV-5297356/relate/chr8-7301025-INV-5297356.filter.poplabels","chr8-7301025-INV-5297356","output/chr8/chr8-7301025-INV-5297356/relate/all.nwk","chr8-7301025-INV-5297356")

fin_tped <- args[1]
fin_hapID <- args[2]
INV_ID <- args[3]
fin_nwk <- args[4]
#fin_timetree <- args[5]
outf_prefix <- args[5]

#orig1_tped <- read.delim("../../PCA_hap/chr11_50136371_INV_206505.filter.noMissing.vcf.hap.tped", header=F, sep='')
#orig_tfam <- read.delim("../../PCA_hap/chr11_50136371_INV_206505.filter.noMissing.vcf.hap.tfam", header=F, sep='')

orig1_tped <- read.delim(fin_tped, header=F, sep='')
orig_tped <-  t(orig1_tped[,-c(1,2,3,4,5)])

orig_tfam <- read.delim(fin_hapID, sep='')

orig_tfam$V7 <- apply(orig_tfam, 1, function(x) paste(unlist(strsplit(as.character(x),"_"))[2:3],collapse="_"))



## create a mapping between hapID and INV status
tmp_hap1 <- rep(paste(orig_tfam$V7, "1", sep='_'), each=2)
tmp_hap2 <- rep(paste(orig_tfam$V7, "2", sep='_'), each=2)

tmp_hap3 <- rep(paste(orig_tfam$sample, "1", sep='_'), each=2)
tmp_hap4 <- rep(paste(orig_tfam$sample, "2", sep='_'), each=2)

num_sample <- dim(orig_tfam)[1]
hapID <- rep("", num_sample)
hapID[ seq(1, num_sample*2, 2)] <- tmp_hap1[ seq(1, num_sample*2, 2)]
hapID[ seq(2, num_sample*2, 2)] <- tmp_hap2[ seq(2, num_sample*2, 2)]

orig_hapID <- rep("", num_sample)
orig_hapID[ seq(1, num_sample*2, 2)] <- tmp_hap3[ seq(1, num_sample*2, 2)]
orig_hapID[ seq(2, num_sample*2, 2)] <- tmp_hap4[ seq(2, num_sample*2, 2)]

index_INVsite <- which(grepl(INV_ID,orig1_tped$V2))
index_INVhap <- which(orig_tped[,index_INVsite] == 1)

mapping_hap_INV <- data.frame(hapID= hapID, INV=FALSE, orig_hapID= orig_hapID)
mapping_hap_INV$INV[index_INVhap] <- TRUE


#orig_hapID <- read.delim("IDs_wPopID.hapID", header=F, sep='')
#orig_hapID <- read.delim(fin_hapID, header=F, sep='')
row.names(orig_tped) <- orig_hapID


#iqtree <- read.iqtree("chr11_50136371_INV_206505.treefile" )
# plot unrooted phylogeny
tree <- read.tree(fin_nwk)
num_trees <- length(tree)

if (length(tree) > 1000){
  tree <- sample(tree, 500)
}

  
num_seqs <- length(tree[[1]]$tip.label)

df_tree <- as.data.frame(as_tibble(tree[[1]]))

Pop <- sapply(strsplit(df_tree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_tree)[1]-num_seqs))

df_tree$Pop <- Pop
#df_iqtree$Pop <- ifelse(df_iqtree$Pop == "CMP", "Chimpanzee", df_iqtree$Pop)

#df_iqtree$branchGroup <- ifelse(df_iqtree$label == "CMP", 2, 1)
#df_iqtree$branchGroup <- as.factor(df_iqtree$branchGroup)
df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]
#df_iqtree[df_iqtree$label == "CMP", "INV"] = FALSE

df_tree$hapID <- mapping_hap_INV[ match(df_tree$label, mapping_hap_INV$orig_hapID), "hapID"]

#iqtree2 <- iqtree
#iqtree@phylo$edge.length.mod <- iqtree@phylo$edge.length
#iqtree@phylo$edge.length[1] = iqtree@phylo$edge.length.mod[1]*0.1

func_fort <- function(tr) tr %>% fortify %>% mutate(colID="b")
trees <- lapply(tree, func_fort)
trees[[length(trees)]]$colID <- "a"

options(scipen = 3)

revts_tree <- function(tr){
  x <- tr$x
  mx <- max(x, na.rm=T)
  tr$x <- x - mx
  tr$branch <- tr$branch -mx
  return(tr)
}


new_trees <- lapply(trees, revts_tree)

upXlim <- min(unlist(lapply(new_trees, function(x) min(x$x)))) * 1.1
lowXlim <- min(unlist(lapply(new_trees, function(x) min(x$x)))) * 0.1
list_breaks <- round(seq(upXlim,0,length.out = 4), digits=-4)

  
ggdensitree(new_trees, alpha=0.3, layout = "rectangular", aes(colour=colID)) %<+% df_tree + 
#  ggtree::theme_tree2() +
  geom_tippoint(aes(color=Pop, shape=INV), size=4.5) + scale_shape_manual(values=c(20,17)) + scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS"), labels=c("Afircan","NativeAmerican","EastAsian","European","SouthAsian")) +  
  geom_tiplab(aes(label=hapID, subset=INV), linetype=NA, size=4, hjust=-0.1) +
  scale_x_continuous(limits = c(upXlim,-lowXlim)) + 
#  scale_x_continuous(limits = c(upXlim,-lowXlim), breaks=list_breaks, labels=-list_breaks) + 
  guides(linetype="none", size="none", color=guide_legend(title="Population",nrow=2), shape=guide_legend(nrow=2,title=paste("Inversion",outf_prefix, sprintf("#Local trees (Relate):%s",num_trees),sep="\n"))) + theme(legend.position="top", legend.text = element_text(size=16), legend.title = element_text(size=16), legend.key.size = unit(1, 'cm'), axis.text.x = element_blank() ) 

ggsave(paste(outf_prefix,".densitree.pdf",sep=''), height = 9, width = 12, limitsize=F,useDingbats=FALSE)


####################

library(TreeDist)
library(circlize)
library(reshape2)
library(ggplot2)
library(dplyr)
require(scales)
library(RColorBrewer)
library(data.table)
library(cowplot)
library(glue)
library(parallel)


ncolors  = 11
get_colors = function(sdf){
#  bot = 0; top = 1
#  breaks = unique( c(quantile(sdf$value, probs = seq(0, 1, by = 1/ncolors), na.rm=T)) )
  breaks <- seq(-0.099, 1.001, by = 0.1)
  labels = seq(length(breaks)-1)
  # corner case of only one %id value
  if( length(breaks) == 1){
    return(factor(rep(1, length(sdf$value))))
  }
  return( cut(sdf$value, breaks=breaks, labels = labels, include.lowest=TRUE)  )
}


diamond <- function( row ){
#  side_length = as.numeric(row["window"])
  side_length = sqrt(2)
  x=as.numeric(row["w"])
  y=as.numeric(row["z"])
  
  base <- matrix(c(1, 0, 0, 1, -1, 0, 0, -1), nrow = 2) * sqrt(2) / 2
  trans <- (base * side_length) + c(x,y)
  df = as.data.frame(t(trans))
  colnames(df) = c("w","z")
  df$discrete = as.numeric(row["discrete"])
  df$group = as.numeric(row["group"])
  df$value = as.numeric(row["value"])
  df
}

# get the lowest 0.1% of the data so we can not plot it
make_k = function(vals){
  comma(vals/1e3)
}

make_hist = function(sdf, listColors){
  scale_x_fun <- function(x) sprintf("%.1f", x)
  
  bot = quantile(sdf$value, probs=0.001, na.rm=T)[[1]]
  count = nrow(sdf)
  extra = ""
  my_scale = comma
#  if(count > 1e5){
#    extra = "\n(thousands)"
#    my_scale = make_k
#  }
#  p = ggplot(data=sdf, aes(value, fill = discrete)) + geom_histogram(bins=50) +
  p = ggplot(data=sdf, aes(value, fill = factor(discrete))) + geom_histogram(bins=100) +
    theme_cowplot() +
    scale_fill_manual(values=listColors) +
#    scale_fill_brewer(palette = "Spectral", direction = -1) + 
    theme(legend.position = "none") + 
    scale_y_continuous(labels=my_scale) +
    scale_x_continuous(labels=scale_x_fun) +
    coord_cartesian(xlim = c(0, 1))+
    xlab("% similarity")+ylab(glue("# tree pairs {extra}"))
  p
}


make_tri = function(ssdf, var_ID="", listColors){
  ssdf$w = (ssdf$first_pos  + ssdf$second_pos) 
  ssdf$z = -ssdf$first_pos + ssdf$second_pos
#  tri_scale =  max(ssdf$first_pos)/max(ssdf$w) 
#  ssdf$window = 1/tri_scale
#  tri_scale =   1/2.8
 # ssdf$window = 2.8
#  window = 2.8
  ssdf$group = seq(nrow(ssdf))
  df.d = rbindlist(apply(ssdf, 1, diamond ))
  numWindows = sqrt(dim(ssdf)[1])
  
  ggplot(df.d)+
#    geom_polygon(aes(x = w*tri_scale, y = z*window , group=group, fill = factor(discrete))) + 
    geom_polygon(aes(x = w, y = z, group=group, fill = factor(discrete))) + 
    theme_cowplot() +
    scale_fill_manual(values=listColors, na.value = "white") + 
#    scale_fill_gradient2(low = "steelblue", mid ="peru", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1)) +
#    scale_fill_brewer(palette = "Spectral", direction = -1) + 
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(0, NA)) +
    xlab("") + ylab("") +
    theme(legend.position = "none", 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()) + 
    ggtitle(paste(var_ID," (",numWindows," windows [Relate])",sep=""))
}


# main plotting code
col_fun<-colorRamp2(c(0,0.5,1), c("darkblue","white", "darkred"))

mycolors_tri <- colorRampPalette(c("gray","peru","darkred"))(12)
names(mycolors_tri) <- seq(12)

#mat_mutualClustering <- MutualClusteringInfo(tree, normalize = T)
#mat_mutualClustering <- as.matrix(1 - wRF.dist(tree, normalize = T))
mat_mutualClustering <- MatchingSplitInfo(tree, normalize = T)
diag(mat_mutualClustering) <- NA
mat_mutualClustering[ lower.tri(mat_mutualClustering)] <- NA
df_mat_mutualClustering <- reshape2::melt(mat_mutualClustering)
df_mat_mutualClustering$discrete <- get_colors(df_mat_mutualClustering)
names(df_mat_mutualClustering) <- c("first_pos","second_pos","value","discrete")

sdf <- df_mat_mutualClustering
p_lone = make_tri(sdf, var_ID=INV_ID, listColors = mycolors_tri)
scale = 1/2.25
p_hist = make_hist(sdf, listColors = mycolors_tri)
p_hist = p_hist + 
  theme(text = element_text(size=14), axis.text=element_text(size=12))

# ranges for inset hist
mmax = max(sdf$first_pos, sdf$second_pos)
build = ggplot_build(p_lone)
yr = build$layout$panel_params[[1]]$y.range
xmin = - 1/10 * mmax
xmax = mmax * 1/2
ymin = yr[2] * 1/1.8
ymax = yr[2] * 2/2
print(paste(INV_ID, xmin, xmax, ymin, ymax))
# combine
plt = p_lone + annotation_custom(ggplotGrob(p_hist),
  xmin = xmin, xmax = xmax, 
  ymin = ymin, ymax = ymax
)

ggsave(plot=plt,file=glue("{outf_prefix}.TreeSimilarity.pdf"),
#       height=12*scale,
       height=6,
       width = 9, useDingbats=FALSE)


#p1<-ggplot(df_mat_mutualClustering, aes(Var1, Var2, fill=value)) + geom_tile() + theme_void() + scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1))

#p1<-ggplot(rotate(df_mat_mutualClustering,225), aes(Var1, Var2, fill=value)) + geom_tile() + theme_void() + scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1))
#print(p1, vp=viewport(angle=90))


#p1 <- Heatmap(mat_mutualClustering, cluster_rows = F, cluster_columns = F, col=col_fun)
