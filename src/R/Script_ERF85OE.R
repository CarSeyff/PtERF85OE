setwd("~/Git/UPSCb/projects/ERF/ERF85-2018")

# to variance stabilise the data
suppressPackageStartupMessages(library(DESeq2))

# to check the VST fit
suppressPackageStartupMessages(library(vsn))

# to plot the heatmap
suppressPackageStartupMessages(library(gplots))

# for pretty colors
suppressPackageStartupMessages(library(RColorBrewer))

# for reading kallisto data
suppressPackageStartupMessages(library(tximport))

# for plotting the mean-variance relationship
suppressPackageStartupMessages(library(vsn))

#' For plotting scatterplots in 3D
suppressPackageStartupMessages(library(scatterplot3d))

#' Define a palette
pal <- brewer.pal(8,"Dark2")

#' Save the default margin parameters
mar=par("mar")

#' source a few helper scripts
source("~/Git/UPSCb/src/R/featureSelection.R")
source("~/Git/UPSCb/src/R/plot.multidensity.R")

counts <- read.csv("C://Users/Carolin/OneDrive - Umeå universitet/1_Projects_Caro/ERFs/ERF85/BioinformaticPart/ResequencedData2018/HTSeq/HTSEq_ERF85r.txt", sep="\t", header = T, row.names = 1)
#samples<-read.csv("samples.txt", sep="\t", header = T)

counts <- as.matrix(counts)
storage.mode(counts) = "integer"
head (counts)
meta<-data.frame(row.names=colnames(counts),condition=c(rep("T89",3), rep("ERF85",3)))

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = data.frame(condition=meta$condition),
  design = ~ condition)

dds <- estimateSizeFactors(dds)
sizes <- sizeFactors(dds)
names(sizes) <- colnames(counts)
library(pander)
pander(sizes)
boxplot(sizes, main="Sequencing libraries size factor")

#Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))
[1] "21.6% percent (8944) of 41335 genes are not expressed"

#Display the per-gene mean expression i.e. the mean raw count of every gene across samples is calculated and displayed on a log10 scale.
plot(density(log10(counts[,1])))
lines(density(log10(counts[,2])),col="#4d4d4d")
lines(density(log10(counts[,3])),col="#01665e")
lines(density(log10(counts[,4])),col="#35978f")
lines(density(log10(counts[,5])),col="#80cdc1")
lines(density(log10(counts[,6])),col="#8c510a")
abline(v=log10(10))

sel <- rowSums(sapply(lapply(
  split.data.frame(t(counts >= 2),conditions)
  ,colSums), ">=", 2)) >= 1

#The cumulative coverage in our samples is relatively good
#sequencing depth that we achieved is ok 

library("DESeq2")
library("pander")
library("vsn")

vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
colnames(vst) <- colnames(counts)
vst <- vst - min(vst)
write.table(vst, file="vst_ERF85.txt", row.names=T, col.names=T, sep="\t", quote=F)
#Validate the VST
meanSdPlot(vst[rowSums(counts)>0,])
#Visualize the corrected mean - sd relationship. 
#First perform a Principal Component Analysis (PCA) 
pc <- prcomp(t(vst))
summary(pc)

percent <- round(summary(pc)$importance[2,]*100)

#Plot the PCA first two dimensions
pal <- brewer.pal(8,"Dark2")
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(meta$condition)],
     main="Principal Component Analysis",sub="variance stabilized counts")


# do boxplot with PC1 and PC2 to see differences between samples
data_for_plot <- cbind(pc$x, meta)
data_for_plot <- as.data.frame(data_for_plot)

boxplot <- ggplot(data_for_plot, aes(x = condition, y= PC1)) #change y for RDA2
z <- boxplot + geom_boxplot() + scale_x_discrete("")
z

library(corrplot)
cor = cor(counts)
corrplot.mixed(cor(counts, method = "pearson"))
col1 <- colorRampPalette(c("navy", "white", "#b2182b"))
corrplot(cor, order="alphabet", method="circle", col=col1(10), title="Correlation TW_NW_ACC_Mock", tl.pos="lt", type="upper", tl.col="black", tl.cex=0.6, tl.srt=45, addCoef.col="black", pch.cex=0.1, insig="blank")

library(ape)
data.dist = as.dist(1 - cor)
out = hclust(data.dist, method="complete")
plot(out)
plot(as.phylo(out), type = "fan")

##Differential expression and Histogram on log2 and p-Values

sum(rowSums(is.na(counts)) > 0)
cdsFull<-newCountDataSet(counts,as.character(meta$condition))
cdsFull <- estimateSizeFactors(cdsFull)
sizeFactors(cdsFull)
cdsFullBlind <- estimateDispersions(cdsFull, method="blind")
cdsFull <- estimateDispersions(cdsFull)
vsdFull <- varianceStabilizingTransformation(cdsFullBlind)
res = nbinomTest(cdsFull, "T89", "ERF85")
write.table(res, file="res_T89_ERF85.txt", row.names=T, col.names=T, sep="\t", quote=F)

plotMA(res)
plot(res$baseMean,res$log2FoldChange,log="x",pch=ifelse(res$padj<.05,19,20),col=ifelse(res$padj<.05,"red","black"))
hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
res <- subset(res, padj < 0.05)
DEGs_T89_ERF85_0.5fold <- subset(res, res[, "log2FoldChange"] < -0.5 | res[, "log2FoldChange"] > 0.5)
write.table(DEGs_T89_ERF85_0.5fold, file="DEGs_T89_ERF85_0.5fold.txt", row.names=T, col.names=T, sep="\t", quote=F)

### get the vst expression values for your DRGs####
vst <- data.frame(vst)
vst$T89 <- (vst$T89.1 + vst$T89.2 + vst$T89.3)/3
vst$ERF85 <- (vst$X250.2 + vst$X250.5 + vst$X250.9)/4

vst_DEGs_ERF85OE_0.5fold  <- merge(DEGs_T89_ERF85_0.5fold, vst, by.x = "id", by.y = "row.names")

library(pheatmap)

heat_DRGs_ERF85OE <- vst_DEGs_ERF85OE_0.5fold[,c("id", "T89", "ERF85")]
row.names(heat_DRGs_ERF85OE) <- heat_DRGs_ERF85OE[,1]
heat_DRGs_ERF85OE <- heat_DRGs_ERF85OE[,-1]
pheatmap(as.matrix(heat_DRGs_ERF85OE))

##############################################################################################
##############################################################################################
##############################################################################################
