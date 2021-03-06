####################################################################
# Author    : M. Dubbelaar
# Date      : 30-sept-2015
# File Name : DESeq2Visualisation.R
# Purpose   : Makes a visualisation of the APP23 data with the 
#             package DESeQ2.
# Used Files: DESeq2Functions.R
#             DESeq2DE.R
#             DESeq2LinearTime.R
#             saveDEData.R
####################################################################
#                     Loading necessary items                      #
####################################################################
setwd("/home/mdubbelaar/APP23/RScripts/DEseq2/")
source("http://bioconductor.org/biocLite.R")
#biocLite("biomaRt")
#biocLite("DESeq2")
#install.packages("RColorBrewer")
#biocLite("gplots")
library("biomaRt")
library("DESeq2")
library("RColorBrewer")
library("gplots")
resultPathway <- "/home/mdubbelaar/Desktop/APP23_results/DEseq2/"
####################################################################
#                            Loading data                          #
####################################################################
directory <-  "/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/perSampleExpression/"
sampleFiles <- grep("sample", list.files(directory), value = T)
targets <- read.table("/home/mdubbelaar/APP23/Targets.csv", sep=",", header = T)
source("Code/DESeq2Functions.R")
####################################################################
#                    Differential Expression                       #
####################################################################
source("Code/DESeq2DE.R")
####################################################################
#                   Visualisation of the data                      #
####################################################################
# The results are saved within the vector res
# res shows the basemean, the log2FoldChange, the lfcSE, the stat, 
# pvalue and the padj
res <- results(dds)
# A MA plot is made, it shows the log2 fold changes (M) against the 
# mean of the normalized counts (A)
pdf(paste(resultPathway, "Plots/MA_plot.pdf", sep="")) 
plotMA(dds, ylim=c(-3,3), main="DESeq2 MA plot")
dev.off()
####################################################################
# the names of the dds data is changed into the name of the samples 
# to make it easier to check the expression of each sample.
colnames(dds) <- targets$Samples

# The regularized log transformation is usefull to check
rld <- rlogTransformation(dds)
# Calculates the variance stabilizing tranformation and transforms 
# the count data.
vst <- varianceStabilizingTransformation(dds)
# The select vector contains the 30 most highly expressed genes.
select <- order(rowMeans(counts(dds, normalized=T)), decreasing = T)[1:30]
# hmcol makes sure that the heatmaps get a nice color.
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

pdf(paste(resultPathway, "Plots/heatmaps_30_mostExpressedGenes.pdf", sep="")) 
# The heatmap shows data of raw counts of the 30 most expressed genes
heatmap.2(counts(dds, normalized=T)[select,], col=hmcol,
          scale = "row", trace = "none", margin=c(8,9), cexRow = 0.8, cexCol = 0.6, main="Raw counts")
# The heatmap contains regularized log transformation data of the 30 most expressed genes
heatmap.2(assay(rld)[select,], col=hmcol,
          scale = "row", trace = "none", margin=c(8,9), cexRow = 0.8, cexCol = 0.6, main="Regularized log transformation")
# The heatmap contains variance stabilizing transformation data of the 30 most expressed genes
heatmap.2(assay(vst)[select,], col=hmcol,
          scale = "row", trace = "none", margin=c(8,9), cexRow = 0.8, cexCol = 0.6, main="Var stabelizing transformation")
dev.off()
####################################################################
# outliners and as input for machine learning techniques.
# The function is used to visualize the create a sample-to-sample matix
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
                                       paste(condition))
# hclust performs a hierarchical cluster analysisi on a set of 
# dissimilatities and methods for analysing it.
# the function on hclust within the sample-to-sample plot
# is to cluster the samples together.
hc <- hclust(distsRL)
pdf(paste(resultPathway, "Plots/Sample-to-sample_distances.pdf", sep="")) 
heatmap.2(mat, Rowv = as.dendrogram(hc),
          symm = T, trace = "none", col= rev(hmcol))
dev.off()
####################################################################
# The plotPCA helps to check for batch effects.
pdf(paste(resultPathway, "Plots/PCA_plot.pdf", sep="")) 
plotPCA(rld)
dev.off()
####################################################################
# The dispersion estimates are plotted in the plot below.
plotDispEsts(dds)
####################################################################
#                      Benjamini-Hochberg Plot                     #
####################################################################
# The Benjamini-Hochberg multiple testing adjustment procedure graphical illustration.
use <- res$baseMean 
resFilt <- res[use & !is.na(res$pvalue),]
orderInPlot <- order(resFilt$pvalue)
showInPlot <- (resFilt$pvalue[orderInPlot]) <= 0.08
alpha <- 0.05

# The black line shows the p-values vs their rank.
# The smallest p-value can be found in the bottom left corner.
# The red line is the slope (alpha (False Discovery Rate)/ n (the number of tests))
pdf(paste(resultPathway, "Plots/Benjamini-Hochberg_plot.pdf", sep="")) 
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)
dev.off()
####################################################################
#                     Calculate linear timepoints                  #
####################################################################
source("Code/DESeq2LinearTime.R")
setwd("/home/mdubbelaar/Desktop/APP23_results/DEseq2/")
source("/home/mdubbelaar/APP23/RScripts/saveDEData.R")
