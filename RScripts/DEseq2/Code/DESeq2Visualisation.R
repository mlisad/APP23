####################################################################
# Author    : M. Dubbelaar
# Date      : 30-sept-2015
# File Name : DESeq2Visualisation.R
# Purpose   : Makes a visualisation of the APP23 data with the 
#             package DESeQ2.
# Used Files: DESeQ2DE.R
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
####################################################################
#                        Necessary functions                       #
####################################################################
normalizeData <- function(ddsHETSeq) {
  # Normalizes the data.
  # The estimateSizeFactors estimates the size of the given dataset.
  # After this the data will be filtered.
  # It retreives the number of columns and checks if the 
  # expression of each genes is more than 10 for at least 2 genes.
  # The found genes will be used, the other genes will be filtered.
  dds <- estimateSizeFactors(ddsHETSeq)
  #dds <- dds[rowSums(counts(dds)) > 0,]
  nc <- counts(dds, normalized=T)
  filter <- rowSums(nc >= 10) >=2
  dds <- dds[filter,]
}

getUniqueGenes <- function(data) {
  # Gets the unique genes from the DE
  # The genes within the data are tested for an adjusted p value of 0.05
  # These genes are saved and will be returned.
  uniqueGenes <- which(data$padj < 0.05)
  uniqueGenes <- data[uniqueGenes,]
  return(uniqueGenes)
}

calculateDds <- function(sampleTable, fileName1, fileName2, fileName3, dataType) {
  # This function calculates the dds.
  # DESeqDataSetFromHTSeqCount is used to store input values, intermediate 
  # calculations and results of the DE analysis.
  # The data will be filtered and normalized with the own made normalizeData function
  # The result from this function will be used within the DESeQ function
  # This function performs different statistical steps
  # The dds data will be used for the three other functions to determine the 
  # the main effect and the interaction effect.
  ddsHETSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~genotype*age)
  dds <- normalizeData(ddsHETSeq)
  dds <- DESeq(dds, betaPrior = F)
  
  calculateMainEffectGenotype(dds, fileName1)
  calculateMainEffectAge(dds, fileName2, dataType)
  calculateInteractionEffect(dds, fileName3, dataType)
}

calculateMainEffectGenotype <- function(dds, fileName){
  # The main effect of the genotype is for every age the same
  # This is the reason why this function is a small one.
  # The genotype of the WT is measured against the APP (HET) mice.
  # The results are transcribed into a txt file.
  mainEffectGenotype <- results(dds, name = "genotype_HET_vs_WT", alpha = 0.05)
  mainEffectGenotype <- getUniqueGenes(mainEffectGenotype)
  write.table(rownames(mainEffectGenotype), fileName, row.names = F, col.names=F, sep=", ", quote = F)
}

calculateMainEffectAge <- function(dds, fileName, dataType) {
  # The function calculateMainEffectAge needs to compare different ages.
  # The comparisons are made among all the ages, young mice (6-8 weeks and 12 months)
  # and old mice (12 months, 18 months and 24 months.)
  # The name of the comparison is used to get the genes and the unique
  # genes within these results are saved into a txt file.
  if (dataType != "old") {
    mainEffectAge1 <- results(dds, name = "age_6_vs_2", alpha = 0.05)
    mainEffectAge1 <- getUniqueGenes(mainEffectAge1)
    if (dataType == "all") {
      mainEffectAge2 <- results(dds, name =  "age_18_vs_2", alpha = 0.05)
      mainEffectAge2 <- getUniqueGenes(mainEffectAge2)
      mainEffectAge3 <- results(dds, name =  "age_24_vs_2",  alpha = 0.05)
      mainEffectAge3 <- getUniqueGenes(mainEffectAge3)
      mainEffectAge <- unique(c(row.names(mainEffectAge1), row.names(mainEffectAge2), row.names(mainEffectAge3)))
      write.table(mainEffectAge, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    } else {
      write.table(rownames(mainEffectAge1), fileName, row.names = F, col.names=F, sep=", ", quote = F)
    }
  } else {
    mainEffectAge1 <- results(dds, name =  "age_18_vs_6", alpha = 0.05)
    mainEffectAge1 <- getUniqueGenes(mainEffectAge1)
    mainEffectAge2 <- results(dds, name =  "age_24_vs_6",  alpha = 0.05)
    mainEffectAge2 <- getUniqueGenes(mainEffectAge2)
    mainEffectAge <- unique(c(row.names(mainEffectAge1), row.names(mainEffectAge2)))
    write.table(mainEffectAge, fileName, row.names = F, col.names=F, sep=", ", quote = F)
  }
}

calculateInteractionEffect <- function(dds, fileName, dataType) {
  # The function calculateInteractionEffect needs to compare different ages.
  # The comparisons are made among all the ages, young mice (6-8 weeks and 12 months)
  # and old mice (12 months, 18 months and 24 months.)
  # This function is larger because of the different results in each age.
  # This data is also transcribed into a file.
  if (dataType != "young") {
    interaction1 <- results(dds, name="genotypeHET.age18", alpha=0.05)
    interaction1 <- getUniqueGenes(interaction1)
    interaction2 <- results(dds, name="genotypeHET.age24", alpha=0.05)
    interaction2 <- getUniqueGenes(interaction2)
    if (dataType == "all") {
      interaction3 <- results(dds, name="genotypeHET.age6", alpha=0.05)
      interaction3 <- getUniqueGenes(interaction1)
      interaction <- unique(c(row.names(interaction1), row.names(interaction2), row.names(interaction3)))
      write.table(interaction, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    } else {
      interaction <- unique(c(row.names(interaction1), row.names(interaction2)))
      write.table(interaction, fileName, row.names = F, col.names=F, sep=", ", quote = F)
    }
  } else {
    interaction <- results(dds, name="genotypeHET.age6", alpha=0.05)
    interaction <- getUniqueGenes(interaction)
    write.table(rownames(interaction), fileName, row.names = F, col.names=F, sep=", ", quote = F)
  }
}
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
source("../saveDEData.R")
