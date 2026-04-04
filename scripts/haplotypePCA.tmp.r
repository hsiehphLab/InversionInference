#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(gmodels)
library(irlba)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(missMDA)

t_col <- function(color, percent = 70, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  lt.col = NULL
  for (i in 1:length(color)){
    ## Get RGB values for named color
    rgb.val <- col2rgb(color[i])
    
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
                 max = 255,
                 alpha = (100 - percent) * 255 / 100,
                 names = name)
    lt.col <- c(lt.col, t.col)
  }
  ## Save the color
  invisible(lt.col)
}


##### Input data


#args = c("../candidates_above2Mb_auto_100kbChunks_chr8only/output/chr8/chr8-7/7-63356276-INV-31891/data/chr7-63356276-INV-31891.filter.noMissing.vcf.hap.tped",  "output/chr7/chr7-63356276-INV-31891/data/IDs_wPopID.hapID", "chr7-63356276-INV-31891", "output/chr7/chr7-63356276-INV-31891/data/chr7-63356276-INV-31891.filter.noMissing.vcf.hap.tfam")


#args = c("output/chr8/chr8-7301025-INV-5297356/data/chr8-7301025-INV-5297356.filter.noMissing.vcf.hap.tped","output/chr8/chr8-7301025-INV-5297356/data/IDs_wPopID.hapID","chr8-7301025-INV-5297356")

args = c("output/chr11/chr11-50136371-INV-206505/data/chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tped", "output/chr11/chr11-50136371-INV-206505/data/IDs_wPopID.hapID","chr11-50136371-INV-206505", "output/chr11/chr11-50136371-INV-206505/data/chr11-50136371-INV-206505.filter.noMissing.vcf.hap.tfam")

#args = c("output/chr3/chr3-195749464-INV-230745/data/chr3-195749464-INV-230745.filter.noMissing.vcf.hap.tped", "output/chr3/chr3-195749464-INV-230745/data/IDs_wPopID.hapID","chr3-195749464-INV-230745", "output/chr3/chr3-195749464-INV-230745/data/chr3-195749464-INV-230745.filter.noMissing.vcf.hap.tfam")

fin_tped = args[1]
fin_hapID = args[2]
INVID = args[3]
fin_tfam = args[4]
outf_prefix = args[5]

#orig1_tped <- read.delim("chr11_50136371_INV_206505.filter.noMissing.vcf.hap.tped", header=F, sep='')
orig1_tped <- read.delim(fin_tped, header=F, sep='', na.strings = ".")
orig_tped <-  t(orig1_tped[,-c(1,2,3,4)])

index_INVsite <- which(grepl(INVID,orig1_tped$V2))
index_INVhap <- which(orig_tped[,index_INVsite] == 1)

#orig_hapID <- read.delim("IDs_wPopID.hapID", header=F, sep='')
orig_hapID <- read.delim(fin_hapID, header=F, sep='')
#row.names(orig_tped) <- orig_hapID[,2]

tped <- orig_tped
hapID <- orig_hapID

hapID$V2 <- factor(hapID$V2)


hapID$Pop <- factor(hapID$V2, levels=c("AFR", "AMR", "EAS", "EUR","SAS"))
mycolors_tmp <- colorRampPalette(c('#a50026','#ffeda0','#92c5de','#313695'))(5)

col <- colorRampPalette(mycolors_tmp)(length(unique(hapID$Pop)))[factor(hapID$Pop)]

hapID$indivPop <- sapply(strsplit(as.character(hapID$V1), split="_"), function(x) x[[2]])

new_col <- t_col(col, percent = 50)
names(new_col) <- as.character(hapID$Pop)
names(col) <- as.character(hapID$Pop)
my_pch <- rep(20, length(new_col))
my_pch[index_INVhap] <- 17
my_pch_size <- rep(2, length(new_col))
my_pch_size[index_INVhap] <- 2.5

checkUniqEntry <- function(column){
  l_entry <- unique(column)
  tmp <- unique(l_entry)
  num_nonNA_entry <- length(tmp[!is.na(tmp)])
  return(num_nonNA_entry)
}

uniquelength <- sapply(as.data.frame(tped), checkUniqEntry)
new_tped <- as.data.frame(subset(tped, select=uniquelength>1))
upper_k = min(dim(new_tped)) - 1

if (any(is.na(new_tped))){
  nb <- estim_ncpPCA(new_tped, ncp.min=0, ncp.max=30)
  new_tped.impute <- imputePCA(new_tped, ncp=nb$ncp)
  tped.pca <- prcomp(new_tped.impute$completeObs)
#  tped.pca$rotation

} else{
  tped.pca <- prcomp_irlba(new_tped, n=upper_k)
}

summary_pca <- summary(tped.pca)
summary_pca$importance    
list_propVar_PC <- apply(data.frame(ID=names(summary_pca$importance[2,]), PC=round(summary_pca$importance[2,] * 100, digits=1)), 1, function(x) sprintf("%s(%s%%)",x[1],x[2]))

eigenvec <- tped.pca$x
project.pca <- eigenvec


# k-means
numK <- min(which(cumsum(summary_pca$importance[2,])>0.95)) # find num of PCs explains >85% variance
mydata <- project.pca[,1:numK]
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

for (i in 2:(numK)) wss[i] <- sum(kmeans(mydata, centers=i,  nstart=25, iter.max=1000)$withinss)

#pdf("kmeans_screeplot.pdf",height = 4, width = 6)
pdf(paste(outf_prefix, "kmeans_screeplot.pdf", sep="_"),height = 4, width = 6, useDingbats=FALSE)
plot(1:(numK), wss, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
dev.off()



#### plot pairwise PCs ####
if(numK > 2){
  df_pairwisePC <- combn(1:3, m=2)
}else{
  df_pairwisePC <- combn(1:numK, m=2)
}

for (ii in 1:dim(df_pairwisePC)[2]){
  i=df_pairwisePC[1,ii]; j=df_pairwisePC[2,ii]; pdfFile=sprintf("%s.PC%ivsPC%i.pdf", outf_prefix, i, j)
  pdf(pdfFile, height = 9, width = 9, useDingbats=FALSE)
  plot(project.pca[,i], project.pca[,j], type="n", main='', adj=0.5, xlab=list_propVar_PC[i], ylab=list_propVar_PC[j], bty='n')
  points(project.pca[,i], project.pca[,j], col=new_col, pch=my_pch, cex=my_pch_size)
  legend("bottomright", bty="n", inset=c(0,1), xpd=TRUE, cex=1.0, ncol=5, title="",
         c("AFR", "AMR", "EAS", "EUR","SAS"), fill=mycolors_tmp)
  dev.off()
}

### output some files for later usage
library(factoextra)
iris_transform = as.data.frame(tped.pca$x[,1:numK])

fviz_nbclust(iris_transform, kmeans, method = 'wss', k.max = numK)
ktmp <- fviz_nbclust(iris_transform, kmeans, k.max=numK, print.summary = T)
k = which.max(ktmp$data$y) # based on the two plots above

kmeans_iris = kmeans(iris_transform, centers = k, nstart = 500)
fviz_cluster(kmeans_iris, data = iris_transform[,1:2])


orig_tfam <- read.delim(fin_tfam, header=F, sep='')

orig_tfam$V7 <- apply(orig_tfam, 1, function(x) paste(unlist(strsplit(as.character(x),"_"))[2:3],collapse="_"))


poplabels <- data.frame(sample=orig_tfam$V1, population="POP_0", group="GRP0", sex=NA)
write.table(poplabels, "../data/poplabels.relate", quote = F, row.names = F, sep='\t')

sample_relate1 <- data.frame(ID_1=0, ID_2=0, missing=0)
sample_relate <- rbind(sample_relate1,	data.frame(ID_1=orig_tfam$V1, ID_2=orig_tfam$V1, missing=0))
write.table(sample_relate, "../data/sample.relate", quote = F, row.names = F, sep='\t')


orig_hapID$V3 <- gsub("hap","", orig_hapID$V1)

dateFile <- data.frame(V1=paste("CMP",orig_hapID$V3,sep=","), V2=-6)
write.table(dateFile, "../data/dateFile.txt", quote = F, col.names=F, row.names = F, sep='\t')

df_hclust <- orig_hapID
df_hclust$cluster <- kmeans_iris$cluster
write.table(df_hclust, "../data/df_hclust.txt", quote = F, col.names=F, row.names = F, sep='\t')
#write.table(df_hclust, "output/chr3/chr3-195749464-INV-230745/data/df_hclust.txt", quote = F, col.names=F, row.names = F, sep='\t')


