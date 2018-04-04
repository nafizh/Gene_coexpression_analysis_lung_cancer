scaled.dat <- datExpr


data <- filt_tdat
nsamp <- dim(data)[2] 
data  <- data[,1:nsamp]

# Plot distribution of expression values before normalization
#png("pre-norm.png", 1200, 1200, res=150)
h <- hist((data[,1] ), plot = FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="expression value", ylab="Proportion of molecules")
for(i in 2:nsamp)
{
  h <- hist((data[,i]), plot = FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}
#devnum <- dev.off()

###############################################################################

# MedianNorm function borrowed from the EBSeq library version 1.1.6
# See http://www.bioconductor.org/packages/devel/bioc/html/EBSeq.html
MedianNorm <- function(data)
{
  geomeans <- exp( rowMeans(log(data)) )
  apply(data, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
}

# Normalize by median
size.factors <- MedianNorm(data.matrix(data))
data.norm <- t(apply(data, 1, function(x){ x / size.factors }))

# Plot distribution of normalized expression values
png("post-norm.png", 1200, 1200, res=150)
h <- hist(log(data.norm[,1]), plot=FALSE)
plot(h$mids, h$density, type="l", col=rainbow(nsamp)[1], main="",
     xlab="Log normalized expression value", ylab="Proportion of molecules")
for(i in 2:nsamp)
{
  h <- hist(log(data.norm[,i]), plot=FALSE)
  lines(h$mids, h$density, col=rainbow(nsamp)[i])
}
devnum <- dev.off()


write.table(data, file="log_mean_normalized", sep="\t")
datExpr0 <- read.table("log_mean_normalized")
datExpr0 <- t(datExpr0)
#######################################################################################
# Load the WGCNA package
library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

gsg = goodSamplesGenes(datExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
# outliers.

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#######################################################################################


