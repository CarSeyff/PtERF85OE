#################
################# Script for Bar graphs in Figure 2

library(RColorBrewer)
library(ggplot2)
library(Formula)

######################### 23-12-2019
######################### general graph on Cell area [µm2]
#########################
Stats_CellArea_CellArea<-read.csv("C:/Users/Carolin/OneDrive - Umeå universitet/1_Projects_Caro/ERFs/ERF85/CellWallThickness/NewCalculation_20191223_ERF85OE_cellarea2.txt", sep="\t")

library("ggpubr")
height <- levels(Stats_CellArea$area)
ggboxplot(Stats_CellArea, x = "genotype", y = "area", 
          color = "genotype", 
          ylab = "Fiber cell diameter", xlab = "genotype")

ggviolin(Stats_CellArea, x = "genotype", y = "area", 
         color = "genotype", 
         ylab = "Fiber cell diameter", xlab = "genotype")

qqnorm(Stats_CellArea$area)
qqline(Stats_CellArea$area, col='red')

######################### 02-01-2020
#########################general graph Lumen to lumen length [µm]
#########################
Stats_CellWall<-read.csv("C:/Users/Carolin/OneDrive - Umeå universitet/1_Projects_Caro/ERFs/ERF85/CellWallThickness/vessel_area.txt", sep="\t")

head(Stats_CellWall)

library("ggpubr")
cellwall <- levels(Stats_CellWall$area)
ggboxplot(Stats_CellWall, x = "genotype", y = "area", 
          color = "genotype", 
          ylab = "Lumen to lumen", xlab = "genotype")

ggviolin(Stats_CellWall, x = "genotype", y = "area", 
         color = "genotype", 
         ylab = "Lumen to lumen", xlab = "genotype")

qqnorm(Stats_CellWall$area)
qqline(Stats_CellWall$area, col='red')

#### Statistically it is best to sum up the Replicate in means and to perform statistics using mean values

WT=subset(Stats_CellWall, Stats_CellWall$genotype=="WT")
[750]
Line2=subset(Stats_CellWall, Stats_CellWall$genotype=="Line2")
Line5=subset(Stats_CellWall, Stats_CellWall$genotype=="Line5")
Line9=subset(Stats_CellWall, Stats_CellWall$genotype=="Line9")
Line12=subset(Stats_CellWall, Stats_CellWall$genotype=="Line12")
Line13=subset(Stats_CellWall, Stats_CellWall$genotype=="Line13")

library(RColorBrewer)
library(ggplot2)
library(Formula)

WT_r1=subset(WT, WT$replicate=="Rep1")
WT_r2=subset(WT, WT$replicate=="Rep2")
WT_r3=subset(WT, WT$replicate=="Rep3")

WT_r1_m <- Rle(WT_r1$length)
WT_r2_m <- Rle(WT_r2$length)
WT_r3_m <- Rle(WT_r3$length)

mean(WT_r1_m)
mean(WT_r2_m)
mean(WT_r3_m)

L2_r1=subset(Line2, Line2$replicate=="Rep1")
L2_r2=subset(Line2, Line2$replicate=="Rep2")
L2_r3=subset(Line2, Line2$replicate=="Rep3")

L2_r1_m <- Rle(L2_r1$length)
L2_r2_m <- Rle(L2_r2$length)
L2_r3_m <- Rle(L2_r3$length)

mean(L2_r1_m)
mean(L2_r2_m)
mean(L2_r3_m)

L5_r1=subset(Line5, Line5$replicate=="Rep1")
L5_r2=subset(Line5, Line5$replicate=="Rep2")
L5_r3=subset(Line5, Line5$replicate=="Rep3")

L5_r1_m <- Rle(L5_r1$length)
L5_r2_m <- Rle(L5_r2$length)
L5_r3_m <- Rle(L5_r3$length)

mean(L5_r1_m)
mean(L5_r2_m)
mean(L5_r3_m)

L9_r1=subset(Line9, Line9$replicate=="Rep1")
L9_r2=subset(Line9, Line9$replicate=="Rep2")
L9_r3=subset(Line9, Line9$replicate=="Rep3")

L9_r1_m <- Rle(L9_r1$length)
L9_r2_m <- Rle(L9_r2$length)
L9_r3_m <- Rle(L9_r3$length)

mean(L9_r1_m)
mean(L9_r2_m)
mean(L9_r3_m)

L12_r1=subset(Line12, Line12$replicate=="Rep1")
L12_r2=subset(Line12, Line12$replicate=="Rep2")
L12_r3=subset(Line12, Line12$replicate=="Rep3")

L12_r1_m <- Rle(L12_r1$length)
L12_r2_m <- Rle(L12_r2$length)
L12_r3_m <- Rle(L12_r3$length)

mean(L12_r1_m)
mean(L12_r2_m)
mean(L12_r3_m)

L13_r1=subset(Line13, Line13$replicate=="Rep1")
L13_r2=subset(Line13, Line13$replicate=="Rep2")
L13_r3=subset(Line13, Line13$replicate=="Rep3")

L13_r1_m <- Rle(L13_r1$length)
L13_r2_m <- Rle(L13_r2$length)
L13_r3_m <- Rle(L13_r3$length)

mean(L13_r1_m)
mean(L13_r2_m)
mean(L13_r3_m)

#save in Excel as text

dat <- read.delim("C:/Users/Carolin/OneDrive - Umeå universitet/1_Projects_Caro/ERFs/ERF85/CellWallThickness/Book3.txt", header=T)

qqnorm(dat$lumenToLumen)
qqline(dat$lumenToLumen, col='red')
shapiro.test(dat$area)

[1]p-value = 0.7921 #cell number
[1]p-value =0.9653 #lumen to lumen
[1]p-value =0.1582 #fiber diameter
[2]p-value =0.03204 #vessel diameter
#From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution
#In other words, we can assume the normality.


# For vessel diameter: not normally distributed, use Wilcoxon test (see further below)

genotype <- levels(dat$genotype)
replicate <- levels(dat$replicate)

#Linear model analysis
dat_lmer <- summary(lmer( cell_number~ genotype - 1 + (1|replicate), data=dat))

#Extract means&standard errors of coefficients
estim <- dat_lmer$coefficients[,1]

#Calculate of degree of freedom              
df <- length(dat_lmer$residuals) - length(estim) -  (sum(as.vector(dat_lmer$ngrps)) - 1 )

#Extract variance-covariance matrix
vcov <- as.matrix(dat_lmer$vcov)

#Define genotypeXtreatment interaction
genotype.treatment <- names(estim)

#Make genotypeXtreatment comparison vector 
genotype.treatmentcomp <- c()								
for (i in 1:(length(genotype.treatment)-1)){							
  for (j in (i+1):length(genotype.treatment)){
    genotype.treatmentcomp <- c(genotype.treatmentcomp, paste(genotype.treatment[i], genotype.treatment[j], sep="-"))
  }
}
genotype.treatmentcomp = gsub("genotype", "", genotype.treatmentcomp)


#Calculate p-values (genotype:treatment vs genotype:treatment)
id.mat <- matrix(0, ncol=1, nrow = length(genotype.treatment) )
p.val <- c()
for ( i in 1:(length(estim)-1)) {
  for (j in (i+1):length(estim)){
    id.mat.x <- id.mat
    id.mat.x[ i, 1 ] <- 1
    id.mat.x[ j, 1 ] <- -1
    stder <- sqrt( t(id.mat.x) %*% vcov %*% id.mat.x )
    t.val <- abs( estim[i]-estim[j]) / stder
    p.val <- c( p.val, 2 * pt( t.val, df, lower.tail=F ) )
  }
}
names(p.val) <- genotype.treatmentcomp

library(multcompView)
p.val.lett = multcompLetters(p.val, compare ="<", threshold = 0.05, Letters=c(letters, LETTERS, "."), reverse=F)
p.val.lett

######################### 02-01-2020
#########################Statistics on vessel diameter [µm2]
#########################

se <- function(x) sqrt(var(x)/length(x))

WT=subset(Stats_CellWall, Stats_CellWall$genotype=="WT")
Line2=subset(Stats_CellWall, Stats_CellWall$genotype=="Line2")
Line5=subset(Stats_CellWall, Stats_CellWall$genotype=="Line5")
Line9=subset(Stats_CellWall, Stats_CellWall$genotype=="Line9")
Line12=subset(Stats_CellWall, Stats_CellWall$genotype=="Line12")
Line13=subset(Stats_CellWall, Stats_CellWall$genotype=="Line13")

WT_m <- Rle(WT$area)
Line2_m <- Rle(Line2$area)
Line5_m <- Rle(Line5$area)
Line9_m <- Rle(Line9$area)
Line12_m <- Rle(Line12$area)
Line13_m <- Rle(Line13$area)

mean(WT_m)
mean(Line2_m)
mean(Line5_m)
mean(Line9_m)
mean(Line12_m)
mean(Line13_m)

se(WT_m)
se(Line2_m)
se(Line5_m)
se(Line9_m)
se(Line12_m)
se(Line13_m)


wilcox.test(WT$area,Line2$area)
[p-value = 0.871]

wilcox.test(WT$area,Line5$area)
[p-value = 0.2464]

wilcox.test(WT$area,Line9$area)
[p-value = 0.004535]

wilcox.test(WT$area,Line12$area)
[p-value = 0.6628]

wilcox.test(WT$area,Line13$area)
[p-value = 0.1612]





