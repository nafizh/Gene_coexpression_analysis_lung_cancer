library(survival)
library(dplyr)
lung.dat <- read.table("genomicMatrix_lung_array")
lung.clin.dat <- read.delim("clinical_data_lung_array")

# For clinical data, get only rows which do not have NA in column "X_EVENT"
lung.no.na.dat <- lung.clin.dat[!is.na(lung.clin.dat$X_EVENT), ]

# Getting the transpose of main lung cancer data
ge <- t(lung.dat)
# Getting a vector of all the id's in the clinical data frame without any 'NA' values
keep <- lung.no.na.dat$sampleID

# getting only the samples(persons) for which we have a value rather than 'NA' values
real.dat <- ge[ge[, 1] %in% keep, ]

# adding the 2 columns from clinical data to gene expression data
keep_again <- real.dat[, 1]
temp_df <- lung.no.na.dat[lung.no.na.dat$sampleID %in% keep_again, ]

# naming the columns into our gene expression data
col_names <- ge[1, ]
colnames(real.dat) <- col_names

dd <- temp_df[, c('X_TIME_TO_EVENT', 'X_EVENT')]
real.dat <- cbind(real.dat, dd)

# renaming the 2 new added columns
colnames(real.dat)[colnames(real.dat) == 'X_TIME_TO_EVENT'] <- 'time'
colnames(real.dat)[colnames(real.dat) == 'X_EVENT'] <- 'event'

real.dat[,1] <- NULL # deleting the sample column
df <- real.dat # the main data frame with the 'time' and 'event' column
indx <- sapply(df, is.character)
df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))


# filtering to get the low values out of the data set
wdf <- df[,1:12042]
tdf <- as.data.frame(t(wdf)) # converted to data frame so that dplyr can be applied
#filt_tdat <- filter(tdf, rowMeans(tdf) > 1) # final matrix dim - 17537 520
filt_tdat <- subset(tdf, rowMeans(tdf) > 1) # final matrix dim - 17537 520


# ########## finally got a dataframe of gene expression data and corresponding pcox data
# kp <- rownames(filt_tdat)
# mm <- m[m[, 1] %in% kp, ]
# rdf <- filt_tdat[rownames(filt_tdat) %in% mm[,1], ] # final dim 17463 520
# colnames(mm)[2] <- c('pcox')

sf <- standardScreeningCensoredTime(df$time, df$event, df[, 1:12042])
tf <- sf[sf$ID %in% rownames(filt_tdat), ]
rdf <- filt_tdat[rownames(filt_tdat) %in% tf[,1], ]
ndf <- t(rdf) # use ndf in following analysis - gene names in columnn names





# read in the R libraries
library(MASS) # standard, no need to install
library(class) # standard, no need to install
library(cluster)
library(impute)# install it for imputing missing value 
library(WGCNA)
source("NetworkFunctions.txt")
options(stringsAsFactors=F) 

ndf <- t(rdf)
no.samples = dim(ndf)[[1]]
dim(ndf)
gc() 


# To choose a cut-off value, we propose to use the Scale-free Topology Criterion (Zhang and
# Horvath 2005). Here the focus is on the linear regression model fitting index
# (denoted below by scale.law.R.2) that quantify the extent of how well a network
# satisfies a scale-free topology.
# The function PickSoftThreshold can help one to estimate the cut-off value
# when using hard thresholding with the step adjacency function.
# The first column (different from the row numbers) lists the soft threshold Power
# The second column reports the resulting scale free topology fitting index R^2 (scale.law.R.2)
# The third column reports the slope of the fitting line.
# The fourth column reports the fitting index for the truncated exponential scale free model.
# Usually we ignore it.
# The remaining columns list the mean, median and maximum connectivity.
# To a soft threshold (power) with the scale free topology criterion:
# aim for reasonably high scale free R^2 (column 2), higher than say .80
# and negative slope (around -1, col 4).
# In practice, we pick the threshold by looking for a "kink" in the
# relationship between R^2 and power, see below.

#SOFT THRESHOLDING
# Now we investigate soft thesholding with the power adjacency function
powers1=c(seq(1,10,by=1),seq(12,20,by=2))
RpowerTable=pickSoftThreshold(ndf, powerVector=powers1)[[2]] 




gc()
cex1=0.7
par(mfrow=c(1,2))
plot(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],xlab="
     Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
text(RpowerTable[,1], -sign(RpowerTable[,3])*RpowerTable[,2],
     labels=powers1,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.95,col="red")
plot(RpowerTable[,1], RpowerTable[,5],xlab="Soft Threshold (power)",ylab="Mean
     Connectivity", type="n")
text(RpowerTable[,1], RpowerTable[,5], labels=powers1, cex=cex1,col="red") 



# We use the following power for the power adjacency function.
beta1=7
Connectivity=softConnectivity(ndf,power=beta1)-1
# Let’s create a scale free topology plot.
# The black curve corresponds to scale free topology and
# the red curve corresponds to truncated scale free topology.
par(mfrow=c(1,1))
scaleFreePlot(Connectivity, main=paste("soft threshold, power=",beta1), truncated=F)



#Module Detection
# An important step in network analysis is module detetion.
# Here we use methods that use clustering in combination with the topological
# overlap matrix.
# This code allows one to restrict the analysis to the most connected genes,
# which may speed up calculations when it comes to module detection.
ConnectivityCut = 3600 # number of most connected genes that will be considered
# Incidentally, in the paper by Mischel et al (2005) we considered all 3600 #genes.
ConnectivityRank = rank(-Connectivity)
restConnectivity = ConnectivityRank <= ConnectivityCut
# thus our module detection uses the following number of genes
sum(restConnectivity)
# Now we define the adjacency matrix for the 3600 most connected genes
ADJ= adjacency(ndf[,restConnectivity],power=beta1)
gc()
# The following code computes the topological overlap matrix based on the
# adjacency matrix.
# TIME: This about a few minutes....
dissTOM=TOMdist(ADJ)
gc() 
# Now we carry out hierarchical clustering with the TOM matrix.
# This takes a couple of minutes.
hierTOM = hclust(as.dist(dissTOM),method="average");
par(mfrow=c(1,1))
plot(hierTOM,labels=F) 



# According to our definition, modules correspond to branches of the tree.
# The question is what height cut-off should be used? This depends on the
# biology. Large heigth values lead to big modules, small values lead to small
# but tight modules.
# In reality, the user should use different thresholds to see how robust the findings are. 
# The function cutreeStatistColor colors each gene by the branches that
# result from choosing a particular height cut-off.
# GREY IS RESERVED to color genes that are not part of any module.
# We only consider modules that contain at least 125 genes.
colorh1= cutreeStaticColor(hierTOM,cutHeight = 0.98, minSize = 50)
# The above should be identical to colorh1=datSummary$color1[restConnectivity]
par(mfrow=c(2,1),mar=c(2,4,1,1))
plot(hierTOM, main="Cluster Dendrogram", labels=F, xlab="", sub="");
plotColorUnderTree(hierTOM,colors=data.frame(module=colorh1))
title("Module (branch) color") 


#Module significance
#Next we define a gene significance variable as minus log10 of the univarite Cox regression pvalue
#for predicting survival on the basis of the gene epxression info
# this defines the gene significance for all genes
GeneSignificanceALL <- -log10(tf$pvalueDeviance)
# gene significance restricted to the most connected genes:
GeneSignificance <- GeneSignificanceALL[restConnectivity]
# The function verboseBarplot creates a bar plot
# that shows whether modules are enriched with essential genes.
# It also reports a Kruskal Wallis P-value.
# The gene significance can be a binary variable or a quantitative variable.
# also plots the 95% confidence interval of the mean
par(mfrow=c(1,1))
verboseBarplot(GeneSignificance,colorh1,main="Module Significance ",
               col=levels(factor(colorh1)) ,xlab="Module" )


# To get a sense of how related the modules are one can summarize each module
# by its first eigengene (referred to as principal components).
# and then correlate these module eigengenes with each other.
datME=moduleEigengenes(ndf[,restConnectivity],colorh1)[[1]]
# We define a dissimilarity measure between the module eigengenes that keeps track of the sign of
#the correlation between the module eigengenes.
dissimME=1-(t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based on the module eigengenes of modules") 


#Compare this output with the following:
signif(cor(datME, use="p"), 2) 




# Now we create scatter plots of the samples (arrays) along the module eigengenes.
datMEordered=datME[,hclustdatME$order]
pairs( datMEordered, upper.panel = panel.smooth, lower.panel = panel.cor ,
       diag.panel=panel.hist ,main="Relation between module eigengenes") 



#To study how connectivity is related to mean gene expression or variance of gene expression
# we create the following plot.
mean1=function(x) mean(x,na.rm=T)
var1=function(x) var(x,na.rm=T)
meanExpr=apply(ndf[,restConnectivity],2,mean1)
varExpr=apply(ndf[,restConnectivity],2,var1)
par(mfrow=c(1,2))
plot(Connectivity[restConnectivity],meanExpr, col=as.character(colorh1),
     main="Mean(Expression) vs K",xlab="Connectivity")
plot (Connectivity[restConnectivity],varExpr, col= as.character(colorh1), main="Var(Expression)
      vs K" ,xlab="Connectivity")


# The following produces heatmap plots for each module.
# Here the rows are genes and the columns are samples.
# Well defined modules results in characteristic band structures since the corresponding genes are
# highly correlated.
par(mfrow=c(2,1), mar=c(1, 2, 4, 1))
ClusterSamples=hclust(dist(ndf[,restConnectivity] ),method="average")
# for the first (turquoise) module we use
which.module="turquoise"
plotMat(t(scale(ndf[ClusterSamples$order,restConnectivity][,colorh1==which.module ])
),nrgcols=30,rlabels=T, clabels=T,rcols=which.module,
main=which.module )
# for the second (blue) module we use
which.module="blue"
plotMat(t(scale(ndf[ClusterSamples$order,restConnectivity][,colorh1==which.module ])
),nrgcols=30,rlabels=T, clabels=T,rcols=which.module,
main=which.module ) 



# Now we extend the color definition to all genes by coloring all non-module
# genes grey.
color1=rep("grey",dim(ndf)[[2]])
color1[restConnectivity]=as.character(colorh1)
# The function intramodularConnectivity computes the whole network connectivity kTotal,
# the within module connectivity (kWithin). kOut=kTotal-kWithin and
# and kDiff=kIn-kOut=2*kIN-kTotal
ConnectivityMeasures=intramodularConnectivity(ADJ,colors=colorh1)
names(ConnectivityMeasures)
#[1] "kTotal" "kWithin" "kOut" "kDiff"
# The following plots show the gene significance vs intramodular connectivity
colorlevels=levels(factor(colorh1))
par(mfrow=c(3,3),mar=c(5, 4, 4, 2) + 0.1)
for (i in c(1:length(colorlevels) ) ) {
    whichmodule=colorlevels[[i]];restrict1=colorh1==whichmodule
    verboseScatterplot(ConnectivityMeasures$kWithin[restrict1],
                       GeneSignificance[restrict1],col=colorh1[restrict1],main= paste("set I,",
                                                                                      whichmodule),ylab="Gene Significance",xlab="Intramodular k")
} 




# The following data frame contains the kME values corresponding to each module.
datKME=signedKME(ndf, datME)
#Note we have an intramodular connectivity measure for each color.
names(datKME)
#[1] "kMEblue" "kMEbrown" "kMEgreen" "kMEgrey" "kMEturquoise"
#[6] "kMEyellow"
# Note that the intramodular connectivity has been computed for each of the 8000 genes.
dim(datKME)
#[1] 8000 6


#Question: how do the kME measure relate to the standard intramodular connectivity? 
whichmodule="turquoise"
restrictGenes= colorh1== whichmodule
par(mfrow=c(1,1))
verboseScatterplot(ConnectivityMeasures$kWithin[ restrictGenes],
                   (datKME$kMEturquoise[restConnectivity][restrictGenes])^beta1 ,xlab="kIN",ylab="kME^power",
                   col=whichmodule,main="Relation between two measures of intramodular k, ") 



#Question: how do the kME measure relate to the standard intramodular connectivity? 
whichmodule="blue"
restrictGenes= colorh1== whichmodule
par(mfrow=c(1,1))
verboseScatterplot(ConnectivityMeasures$kWithin[ restrictGenes],
                   (datKME$kMEblue[restConnectivity][restrictGenes])^beta1 ,xlab="kIN",ylab="kME^power",
                   col=whichmodule,main="Relation between two measures of intramodular k, ") 

#Question: how do the kME measure relate to the standard intramodular connectivity? 
whichmodule="yellow"
restrictGenes= colorh1== whichmodule
par(mfrow=c(1,1))
verboseScatterplot(ConnectivityMeasures$kWithin[ restrictGenes],
                   (datKME$kMEyellow[restConnectivity][restrictGenes])^beta1 ,xlab="kIN",ylab="kME^power",
                   col=whichmodule,main="Relation between two measures of intramodular k, ") 


#Question: find genes with high gene significance (Cox-pvalue smaller than 0.05) and high
#intramodular connectivity in the blue module. 

FilterGenes <- GeneSignificanceALL > -log10(0.05) & abs(datKME$kMEturquoise) > .85
table(FilterGenes)
gene_tor_sig <- (rdf[FilterGenes,] )

FilterGenes <-  abs(datKME$kMEturquoise) > 0.8
table(FilterGenes)
gene_tor <- (rdf[FilterGenes,] )
final_genes <- rownames(gene_tor)

FilterGenes <-  abs(datKME$kMEbrown) > 0.8
table(FilterGenes)
gene_brown_small <- (rdf[FilterGenes,] )

FilterGenes <-  abs(datKME$kMEblue) > 0.8
table(FilterGenes)
gene_blue <- (rdf[FilterGenes,] )
final_genes <- rownames(gene_blue)

FilterGenes <-  abs(datKME$kMEgreen) > 0.8
table(FilterGenes)
gene_green <- (rdf[FilterGenes,] )

FilterGenes <-  abs(datKME$kMEyellow) > 0.8
table(FilterGenes)
gene_yellow_small <- (rdf[FilterGenes,] )


#array yellow genes were in rna blue 
#array tor genes were in rna tor