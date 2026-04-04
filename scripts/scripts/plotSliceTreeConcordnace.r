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

mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)
cols <- c("steelblue", mycolors_tmp, "black")
#names(cols) <- c("Chimpanzee","AFR","AMR","EAS","EUR","SAS","internal")
names(cols) <- c("b","AFR","AMR","EAS","EUR","SAS","a")
#cols <- mycolors_tmp
#names(cols) <- c("AFR","AMR","EAS","EUR","SAS")



#args = c("output/chrX/chrX-154339398-INV-53354/relate/chrX-154339398-INV-53354.filter.noMissing.anc.haps", "output/chrX/chrX-154339398-INV-53354/data/poplabels.relate", "chrX-154339398-INV-53354", "output/chrX/chrX-154339398-INV-53354/iqtree_sliceFASTA/", "output/chrX/chrX-154339398-INV-53354/iqtree/chrX-154339398-INV-53354.timetree.nex", "output/chrX/chrX-154339398-INV-53354/iqtree_sliceFASTA/all.minMutHomoplasy.dat", "input/refSeq.hg38.wX.uniq.bed", 10000, "tmp", "output/chrX/chrX-154339398-INV-53354/data/locus.bed")

#args = c("output/chr11/chr11-50136371-INV-206505/relate/chr11-50136371-INV-206505.filter.noMissing.anc.haps", "input/poplabels.relate", "chr11-50136371-INV-206505", "output/chr11/chr11-50136371-INV-206505/iqtree_sliceFASTA/", "output/chr11/chr11-50136371-INV-206505/iqtree/chr11-50136371-INV-206505.timetree.nex", "output/chr11/chr11-50136371-INV-206505/iqtree_sliceFASTA/all.minMutHomoplasy.dat", "input/refSeq.hg38.uniq.bed", 20000, "tmp")


#args = c("output/chr8/chr8-7301025-INV-5297356_test1/relate/chr8-7301025-INV-5297356_test1.filter.noMissing.anc.haps", "output/chr8/chr8-7301025-INV-5297356_test1/data/poplabels.relate", "chr8-7301025-INV-5297356_test1", "output/chr8/chr8-7301025-INV-5297356_test1/iqtree_sliceFASTA/", "output/chr8/chr8-7301025-INV-5297356_test1/iqtree/chr8-7301025-INV-5297356_test1.timetree.nex", "output/chr8/chr8-7301025-INV-5297356_test1/iqtree_sliceFASTA/all.minMutHomoplasy.dat", "input/refSeq.hg38.uniq.bed", 20000, "tmp", "output/chr8/chr8-7301025-INV-5297356_test1/data/locus.bed")


#args = c("output/chr17/chr17-45585160-INV-706887/relate/chr17-45585160-INV-706887.filter.haps", "output/chr17/chr17-45585160-INV-706887/relate/chr17-45585160-INV-706887.filter.poplabels", "chr17-45585160-INV-706887", "output/chr17/chr17-45585160-INV-706887/iqtree_sliceFASTA/all.sliceTrees", "tmp")

#args = c("output/chr1/chr1-13084312-INV-62181/relate/chr1-13084312-INV-62181.filter.haps", "output/chr1/chr1-13084312-INV-62181/relate/chr1-13084312-INV-62181.filter.poplabels", "chr1-13084312-INV-62181", "output/chr1/chr1-13084312-INV-62181/iqtree_sliceFASTA/all.sliceTrees", "tmp")


#args = c("output/chr2/chr2-91832041-INV-180624/relate/chr2-91832041-INV-180624.filter.noMissing.anc.haps", "output/chr2/chr2-91832041-INV-180624/data/poplabels.relate", "chr2-91832041-INV-180624", "output/chr2/chr2-91832041-INV-180624/iqtree_sliceFASTA/", "output/chr2/chr2-91832041-INV-180624/iqtree/chr2-91832041-INV-180624.timetree.nex", "output/chr2/chr2-91832041-INV-180624/iqtree_sliceFASTA/all.minMutHomoplasy.dat", "input/refSeq.hg38.wX.uniq.bed", 2000, "tmp", "output/chr2/chr2-91832041-INV-180624/data/locus.bed")

fin_tped <- args[1]
fin_hapID <- args[2]
INV_ID <- args[3]
path_nwk <- args[4]
fin_timetree <- args[5]
fin_homoplasy = args[6]
refseq = args[7]
winSize = as.numeric(args[8])
outf_prefix <- args[9]
locusBED <- args[10]

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


mytimetree2 <- list()
mytimetree2[[INV_ID]] <- treeio::drop.tip(read.nexus(fin_timetree), "CMP")

list_trees <- list.files(path_nwk, full.names = T, pattern = "timetree.nex")

list_trees_start <- unlist(lapply(strsplit(basename(list_trees), split="\\."), function(x) strsplit(x[1],split="-")[[1]][1]))
 
df_list_trees <- data.frame(start=as.numeric(list_trees_start), path=list_trees)
df_list_trees <- df_list_trees[ order(df_list_trees$start),]

mytimetree <- list()
for(f in df_list_trees$path){
  tmpID <- unlist(strsplit(basename(f), split = "\\."))[1]
  mytimetree[[tmpID]] <- treeio::drop.tip(read.nexus(f), "CMP")
}

if (length(mytimetree) > 1000){
  mytimetree <- sample(mytimetree, 500)
}

num_trees <- length(mytimetree)

mytimetree2 <- c(mytimetree2, mytimetree)
class(mytimetree2) <- "multiPhylo" 

sumBranchLength <- function(x, mPhyloObj){
  yrs_gen <- 29
  return(sum(mPhyloObj[[x]]$edge.length*1e6)/yrs_gen)
} 

avail_treeID <- names(mytimetree2)
df_homoplasy <- read.delim(fin_homoplasy, header=F, sep = "")
df_homoplasy <- df_homoplasy[ df_homoplasy$V1 %in% avail_treeID,]

df_homoplasy$totalBranchLength <- apply(df_homoplasy, 1, function(x) sumBranchLength(x[1], mytimetree2))
df_homoplasy$invPerGeneration <- df_homoplasy$V2/df_homoplasy$totalBranchLength


write.table(df_homoplasy, paste(outf_prefix,"df.homoplasySliceTrees.dat", sep = "."))

num_seqs <- length(mytimetree2[[1]]$tip.label)

df_tree <- as.data.frame(as_tibble(mytimetree2[[1]]))

Pop <- sapply(strsplit(df_tree[1:num_seqs,"label"], split = "_"), function(x) x[[1]])
Pop <- c(Pop, rep("internal",dim(df_tree)[1]-num_seqs))

df_tree$Pop <- Pop

df_tree$INV <- mapping_hap_INV[ match(df_tree$label,mapping_hap_INV$orig_hapID), "INV"]

df_tree$hapID <- mapping_hap_INV[ match(df_tree$label, mapping_hap_INV$orig_hapID), "hapID"]


func_fort <- function(tr) tr %>% fortify %>% mutate(colID="b")
trees <- lapply(mytimetree2, func_fort)
trees[[INV_ID]]$colID <- "a"

options(scipen = 3)

revts_tree <- function(tr){
  x <- tr$x
  mx <- max(x, na.rm=T)
  tr$x <- x - mx
  tr$branch <- tr$branch -mx
  return(tr)
}


new_trees <- lapply(trees, revts_tree)

p_densi <- ggdensitree(new_trees, alpha=0.3, tip.order = 1, layout = "rectangular", aes(colour=colID, size=colID)) %<+% df_tree +  scale_size_manual( values = c(2,1) ) +
  ggtree::theme_tree2() +
  geom_tippoint(aes(shape=INV), color="gray10", size=3) + scale_shape_manual(values=c(NA,17), labels=c("","Inversion")) +
  geom_tiplab(aes(color=Pop, label=hapID, subset=INV), key_glyph = "rect", linetype=NA, size=4, hjust=-0.1) +
  guides(linetype="none", size="none", 
         color=guide_legend(title="Population", nrow=2, override.aes = list(size=0.001)), 
         shape=guide_legend(nrow=1, title=paste(outf_prefix, sprintf("#Local trees [slices]:%s",num_trees),sep="\n"))) + 
  scale_color_manual(values=cols, breaks=c("AFR","AMR","EAS","EUR","SAS"), labels=c("Afircan","NativeAmerican","EastAsian","European","SouthAsian")) +
  xlab("Million years") +
  theme(legend.position="top", legend.text = element_text(size=14), legend.title = element_text(size=14), legend.key.size = unit(0.4, 'cm'), axis.text.x = element_text(size=14), axis.title = element_text(size=16)) 


build_densi = ggplot_build(p_densi)
xr = build_densi$layout$panel_params[[1]]$x.range
yr = build_densi$layout$panel_params[[1]]$y.range
xmin = 1/10 * (xr[2] - xr[1])
ymin = 1/10 * (yr[2] - yr[1])

list_breaks <- build_densi$layout$panel_params[[1]]$x.sec$breaks

# basic stats
pointEst_numINVs <- df_homoplasy[ df_homoplasy$V1 == INV_ID, "V2"]
pointEst_invPerGen <- df_homoplasy[ df_homoplasy$V1 == INV_ID, "invPerGeneration"]
mean_numINVs <- mean(df_homoplasy$V2)
CI_numINVs <- quantile(df_homoplasy$V2, probs = c(0.025,0.975))
#sd_numINVs <- sd(df_homoplasy$V2)
#CI_numINVs <- c(mean_numINVs - 1.96*sd_numINVs, mean_numINVs + 1.96*sd_numINVs)
mean_invPerGen <- mean(df_homoplasy$invPerGeneration)
CI_invPerGen <- quantile(df_homoplasy$invPerGeneration, probs = c(0.025,0.975)) 
#sd_invPerGen <- sd(df_homoplasy$invPerGeneration)
#CI_invPerGen <- c(mean_invPerGen - 1.96*sd_invPerGen, mean_invPerGen + 1.96*sd_invPerGen)

annoText1 <- sprintf("Estimated #inversion events:%s (95%% C.I.: [%.2f, %.2f])", pointEst_numINVs, CI_numINVs[1], CI_numINVs[2])
annoText2 <- sprintf("Estimated rate: %.2e (95%% C.I.: [%.2e, %.2e]) inversions/generation", pointEst_invPerGen, CI_invPerGen[1], CI_invPerGen[2])

annoText <- paste(annoText1, annoText2, sep="\n")
annoText.x = xr[1] + xmin * 3.5
annoText.y1 = yr[2] - ymin * 0.5
annoText.y2 = yr[2] - ymin*1

p_densi + annotate("text", label=annoText1, x=annoText.x, y=annoText.y1, size=5) + 
  annotate("text", label=annoText2, x=annoText.x, y=annoText.y2, size=5) + 
  coord_cartesian(xlim=c(NA, xr[2]+xmin)) + scale_x_continuous(breaks=list_breaks,labels=-list_breaks)

ggsave(paste(outf_prefix,".densitree.sliceTrees.pdf",sep=''), height = 9, width = 12, limitsize=F, useDingbats=FALSE)


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
library(ggrepel)

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
#  trans <- base + c(x,y)
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


make_tri = function(ssdf, var_ID="", listColors, y_low = 0, numWindows){
  ssdf$w = (ssdf$first_pos  + ssdf$second_pos) 
  ssdf$z = -ssdf$first_pos + ssdf$second_pos
#  tri_scale =  max(ssdf$first_pos)/max(ssdf$w) 
#  ssdf$window = 1/tri_scale
  ssdf$group = seq(nrow(ssdf))
  df.d = rbindlist(apply(ssdf, 1, diamond ))
  df.d$z <- ifelse(df.d$z < 1, 1, df.d$z)
  ggplot(df.d)+
#    geom_polygon(aes(x = w*tri_scale, y = z*window , group=group, fill = factor(discrete))) + 
    geom_polygon(aes(x = w, y = z, group=group, fill = factor(discrete))) + 
    theme_cowplot() +
    scale_fill_manual(values=listColors, na.value = "grey") + 
#    scale_fill_gradient2(low = "steelblue", mid ="peru", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1)) +
#    scale_fill_brewer(palette = "Spectral", direction = -1) + 
    scale_x_continuous(limits = c(0, NA)) +
    scale_y_continuous(limits = c(y_low, NA)) +
    xlab("") + ylab("") +
    theme(legend.position = "none", 
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank()
    ) + 
    ggtitle(paste(var_ID, " (", numWindows, " windows [slices])", sep=""))
}


# prepare data for plotting genes
INV_coord <- read.delim(locusBED, header=F, sep='')
INV_start <- as.numeric(INV_coord$V2)
INV_end <- as.numeric(INV_coord$V3)
df_refseq <- read.delim(refseq, header = F, sep="")
df_refseq_locus <- df_refseq[ df_refseq$V1 %in% INV_coord$V1 & df_refseq$V2>=INV_start & df_refseq$V3<=INV_end, ] 

window_starts <- seq(INV_start, INV_end, winSize)

if(dim(df_refseq_locus)[1] != 0){
  
  refseq_start_x <- apply(df_refseq_locus, 1, function(x) findInterval(x[2], window_starts))
  refseq_end_x <- apply(df_refseq_locus, 1, function(x) findInterval(x[3], window_starts))
  
  df_refseq_locus_tr <- data.frame(start=refseq_start_x,end=refseq_end_x)
  
  yLow = -(dim(df_refseq_locus_tr)[1] + 1)
  df_refseq_locus_tr$start_w <- 2 * df_refseq_locus_tr$start - 1
  df_refseq_locus_tr$end_w <- 2 * df_refseq_locus_tr$end + 1
  df_refseq_locus_tr$start_z <- seq(-1, yLow+1, -1) 
  df_refseq_locus_tr$end_z <- seq(-1, yLow+1, -1)
  df_refseq_locus_tr$gene <- df_refseq_locus$V4
  
}else{
  df_refseq_locus_tr = data.frame(start=numeric(), end=numeric(), start_w=numeric(), end_w=numeric(), start_z=numeric(), end_z=numeric(), gene=character())
  yLow = 0
}


# main plotting code
col_fun<-colorRamp2(c(0,0.5,1), c("darkblue","white", "darkred"))

mycolors_tri <- c("white", colorRampPalette(c("gray","peru","darkred"))(12))
names(mycolors_tri) <- c(0,seq(12))

#mat_mutualClustering <- MutualClusteringInfo(tree, normalize = T)
#mat_mutualClustering <- as.matrix(1 - wRF.dist(tree, normalize = T))
mat_mutualClustering <- as.matrix(MatchingSplitInfo(mytimetree, normalize = T))

mat_mutualClustering[ is.nan(mat_mutualClustering)] <- NA

mat_idxLookup <- as.data.frame(do.call(rbind, lapply(names(mytimetree), function(x) unlist(strsplit(x,"-")))), stringsAsFactors = FALSE)


mat_idxLookup <- transform(mat_idxLookup, V1 = as.numeric(V1))
mat_idxLookup <- transform(mat_idxLookup, V2 = as.numeric(V2))
mat_idxLookup$scaledV1 <- mat_idxLookup$V1
mat_idxLookup$scaledV2 <- mat_idxLookup$V2


# custom operator for "not in"
'%!in%' <- function(x,y)!('%in%'(x,y))

idx_NA_ScaledMatrix <- which(window_starts %!in% mat_idxLookup[,3])
idx_avail_scaledMatrix <- which(window_starts %in% mat_idxLookup[,3])

idx_origM_scaledMatrix <- rep(NA,length(window_starts))
for (ii in 1:length(idx_avail_scaledMatrix)){
  idx_origM_scaledMatrix[ idx_avail_scaledMatrix[ii]] = ii
}

df_mapping_scaled2origM <- data.frame(V1=1:length(window_starts), V2=idx_origM_scaledMatrix)

mat_mutualClustering_rescaled <- NULL
for (i in 1:length(window_starts)){
  if (i %in% idx_NA_ScaledMatrix){
    mat_mutualClustering_rescaled <- rbind(mat_mutualClustering_rescaled, rep(NA,length(window_starts)))
  }else{
    mat_mutualClustering_rescaled <- rbind(mat_mutualClustering_rescaled, mat_mutualClustering[ df_mapping_scaled2origM[i==df_mapping_scaled2origM$V1,"V2"], idx_origM_scaledMatrix])
  }
}

mat_mutualClustering_rescaled <- as.matrix(mat_mutualClustering_rescaled)
colnames(mat_mutualClustering_rescaled) <- NULL

diag(mat_mutualClustering_rescaled) <- NA
mat_mutualClustering_rescaled[ lower.tri(mat_mutualClustering_rescaled)] <- NA


df_mat_mutualClustering <- reshape2::melt(mat_mutualClustering_rescaled)
df_mat_mutualClustering$discrete <- get_colors(df_mat_mutualClustering)
names(df_mat_mutualClustering) <- c("first_pos","second_pos","value","discrete")

sdf <- df_mat_mutualClustering

options(ggrepel.max.overlaps = Inf)

nWindows = dim(mat_mutualClustering)[1]

if(dim(df_refseq_locus_tr)[1] != 0){
  p_lone = make_tri(sdf, var_ID=INV_ID, listColors = mycolors_tri, y_low = yLow, numWindows = nWindows) + 
    geom_segment(data=df_refseq_locus_tr, aes(x=start_w, y=start_z, xend=end_w, yend=end_z, colour="steelblue"), size=1.5) + 
    scale_color_manual(values="steelblue") + 
    geom_text_repel(data=df_refseq_locus_tr, aes(x=end_w,y=end_z, label=gene))
}else{
  p_lone = make_tri(sdf, var_ID=INV_ID, listColors = mycolors_tri, y_low = yLow, numWindows = nWindows) 
}

#scale = 1/2.25
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


ggsave(plot=plt,file=glue("{outf_prefix}.TreeSimilarity.sliceTrees.pdf"),
      height = 6, width = 9, useDingbats=FALSE)


#p1<-ggplot(df_mat_mutualClustering, aes(Var1, Var2, fill=value)) + geom_tile() + theme_void() + scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1))

#p1<-ggplot(rotate(df_mat_mutualClustering,225), aes(Var1, Var2, fill=value)) + geom_tile() + theme_void() + scale_fill_gradient2(low = "darkblue", mid = "gray", high = "darkred", na.value = "white", midpoint = 0.5, limits=c(-0.1,1.1))
#print(p1, vp=viewport(angle=90))


#p1 <- Heatmap(mat_mutualClustering, cluster_rows = F, cluster_columns = F, col=col_fun)
