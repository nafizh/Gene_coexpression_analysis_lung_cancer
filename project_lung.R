library(survival)
library(dplyr)
lung.dat <- read.table("genomicMatrix_lung_smaller")
lung.clin.dat <- read.delim("clinical_data_lung_smaller")

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
indx <- sapply(df, is.factor)
df[indx] <- lapply(df[indx], function(x) as.numeric(as.character(x)))

datExpr <- df[, 1:20530] # the main data frame without the 'time' and 'event' column


#######################################################################################
# Cleaning and removing genes with 0 variance by WCGNA
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
    datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Next we cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious
# outliers.

sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Checking 2 samples if they are on the 45 degree line
tdat <- t(datExpr)
qqplot(tdat[,1], tdat[,2], col=3)

# MA plot
mm = (tdat[,1]) - (tdat[,2])
aa = (tdat[,1]) + (tdat[,2])
plot(aa,mm,col=2)


# remove low counts - finally kind of fixed the 2 bell curves in the same plot
filt_tdat = as.data.frame(tdat) # converted to data frame so that dplyr can be applied
filt_tdat = filter(filt_tdat, rowMeans(filt_tdat) > 1) # final matrix dim - 17537 520

boxplot(as.matrix((filt_tdat)),col=2)

ematrix = tdat[rowMeans(tdat) > 16, ] # max(rowMeans(tdat)) = 16
colramp = colorRampPalette(c(3,"white",2))(9)
heatmap(ematrix,col=colramp,Rowv=NA,Colv=NA)

heatmap.2(ematrix,col=colramp,Rowv=NA,Colv=NA,
          dendrogram="none", scale="row",trace="none")
#######################################################################################
vec <- rep(NA, 20530)
for(i in 1:ncol(df)){
  cf <- coxph(Surv(time, event) ~ (as.numeric(unlist(df[colnames(df)[i]]))), data = df)
  summcf <- summary(cf)
  vec[i] <- summcf$logtest['pvalue']
} 


m['pcox'] <- vec

vec_new <- vector("list", 20532) # creating an empty list

# Getting the coxph values for all genes
indx1 <- sapply(df, is.double)
vec_new <- lapply(df[indx1], function(x) 
                      summary(coxph(Surv(time, event) ~ (as.numeric(unlist(x))), data = df))$logtest['pvalue']) 


m <- as.data.frame(matrix(0, ncol = 1, nrow = 20530))
colnames(m) <- c('pcox')
i <- c(1:20530)
m$pcox <- lapply(i, function(x) vec_new[[x]]['pvalue'])
mm <- as.data.frame(m)
#####################################################################################


