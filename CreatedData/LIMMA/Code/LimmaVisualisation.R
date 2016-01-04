####################################################################
# Author    : M. Dubbelaar
# Date      : 25-sept-2015
# File Name : LimmaVisualisation.R
# Purpose   : To visualise the APP data with the package limma.
# Used Files: EdgeRFunctions.R
#             loadingEnsemblData.R
#             LimmaLinearTime.R
#             LimmaDispersions.R
#             plotColors.R
#             saveDEData.R
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/APP23/CreatedData/LIMMA/")
#install.packages("gplots")
#install.packages("statmod")
#biocLite("biomaRt")
#biocLite("limma")
#biocLite("edgeR")
library(biomaRt)
library(gplots)
library(limma)
library(edgeR)
library(statmod)
resultPathway <- "/home/mdubbelaar/Desktop/APP23_results/LIMMA/"
####################################################################
#                        Necessary Functions                       #
####################################################################
calculateEffectsLimma <- function(dge, design, namefile1, namefile2, namefile3, pathway) {
  # This function calculates the main effect of genotype, the main effect of age
  # and the interaction between age and genotype. The function voom transforms the 
  # count data to logCPM (log2-counts per million), estimates the main-variance 
  # relationship. The data can be used for linear modelling after these steps.
  # The fit (calculated with the use of the eBayes function) can be used to
  # tests the likelihood of the coefficients. The genes will be filtered 
  # only genes with an adjusted value below the 0.05 will be saved into the
  # result vector. The saved genes will be used in the heatmap for visualisation
  # and the gene names will be saved within a txt file.
  
  # Voom contains several build-in procedures to estimate the main-variance relationship.
  v <- voom(dge, design)
  # lmfit fits the linear model.
  fit <- lmFit(v, design)
  # The fit is made with the use of the eBayes function (note not ebayes this causes errors).
  fit2 <- eBayes(lmFit(v, design))
  
  # The main effects and the interaction effect are defined.
  mainGenotype <- topTable(fit2, coef=2, n=Inf, sort.by="none")
  mainTime <- topTable(fit2, coef = 3, n=Inf, sort.by="none")
  interactionGenotypeTime <- topTable(fit2, coef = 4,  n=Inf, sort.by="none")
  # and the found genes are filtered for an adjusted p value below the 0.05
  mainGenotype.result <- filterGenesWithLimma(mainGenotype)
  mainTime.result <- filterGenesWithLimma(mainTime)
  interactionGenotypeTime.result <- filterGenesWithLimma(interactionGenotypeTime)
  
  # The last step is to plot the genes with their information in a heatmap.
  # The if - else is necessary because of a lack op interaction genes within one
  # of the age groups.
  if (length(rownames(interactionGenotypeTime.result)) == 0) {
    pdf(pathway)
    heatmap.2(M2[rownames(mainTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
    heatmap.2(M2[rownames(mainGenotype.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
    dev.off()  
  } else {
    pdf(pathway)
    heatmap.2(M2[rownames(mainTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
    heatmap.2(M2[rownames(mainGenotype.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
    heatmap.2(M2[rownames(interactionGenotypeTime.result),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="Interaction effect")
    dev.off()
  }
  # The unique genes are saved within a table so these gene names can be used for futher analysis.
  write.table(rownames(mainGenotype.result), namefile1, row.names = F, col.names=F, sep=", ", quote = F)
  write.table(rownames(mainTime.result), namefile2, row.names = F, col.names=F, sep=", ", quote = F)
  write.table(rownames(interactionGenotypeTime.result), namefile3, row.names = F, col.names=F, sep=", ", quote = F)
}

# This function is made to calculate the unique genes for the differential expressions.
# The genes with a adjusted p value < 0.05 are returned.
filterGenesWithLimma <- function(tabel) {
  uniqueGenes <- which(tabel$adj.P.Val < 0.05)
  uniqueGenes <- tabel[uniqueGenes,]
  return(uniqueGenes)
}
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("../EdgeR/Code/EdgeRFunctions.R")
M1 <- getData("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
targets <- getTarget("/home/mdubbelaar/APP23/Targets.csv")
source("../plotColors.R")

# A dge list is creates, this list is used to estimate the GLM 
# common, trended and tagwise dispersion. This object holds the 
# dataset that needs to be analysed by EgdeR and the calculations 
# which where performed on the dataset.
dge <- DGEList(counts=M1)
M2 <- cpm(dge, log=TRUE)
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("../loadingEnsemblData.R")
####################################################################
#                          Filtering Data                          #
####################################################################
# The data is filtered with a cut-off of 1 count per million (cpm).
# The data with 2 replicated samples will be saved into the vector.
isExpr <- rowSums(cpm(dge)>1) >= 2
dge <- dge[isExpr, keep.lib.sizes=F]
# These expressed genes will be saved within the dge and BioM dataset.
BioM <- BioM[isExpr,]
# The raw library sizes are calculated with the use of calNormFactors
# This is a step before calculating the estimates.
dge <- calcNormFactors(dge)
####################################################################
#                             MDS Plot                             #
####################################################################
# The MDS plot is made and saved within a file with the use of pdf()
# dev.off() is necessary to close the connection to the document.
pdf(paste(resultPathway, "Plots/MDS_plot.pdf", sep=""))
plotMDS(dge, labels = targets$Conditie, col=as.numeric(targets$Conditie))
legend("topright", legend = unique(targets$Conditie), col=unique(targets$Conditie), pch=15, cex=.6)
dev.off()
####################################################################
#                              Colors                              #
####################################################################
source("../plotColors.R")
####################################################################
#                      Differential Expression                     #
####################################################################
source("Code/LimmaLinearTime.R")
source("Code/LimmaDispersions.R")
setwd("/home/mdubbelaar/Desktop/APP23_results/LIMMA/")
source("/home/mdubbelaar/APP23/CreatedData/saveDEData.R")
setwd("/home/mdubbelaar/APP23/CreatedData/LIMMA/")
####################################################################
#                            Heatmaps                              #
####################################################################
# Checks the two toptable results with the most unique genes.
# These heatmaps are saved within the heatmap.pdf.
pdf(paste(resultPathway, "Plots/heatmaps_contains_most_unique.pdf", sep=""))
heatmap.2(M2[match(rownames(toptable7.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M HET-18M HET")
heatmap.2(M2[match(rownames(toptable11.results), rownames(M2)),], ColSideColors = col_cell_age, cexRow = 0.01, trace = "none", scale = "row", main="24M HET- 6-8W HET")
dev.off()
