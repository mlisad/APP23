####################################################################
# Author    : M. Dubbelaar
# Date      : 8-dec-2015
# File Name : EdgeRFunctions.R
# Purpose   : This file contains all of the functions that are used
#             in the EdgeRVisualisation file.
####################################################################
#                     Get all of the probe names                   #
####################################################################
setRowname <- function(data) {
  # This function checks for the column with the name "probe"
  # This column will be used as the rownames and the column will be
  # put to null (it will remove the colum).
  colNumber <- grep("probe", colnames(data))
  row.names(data) <- data[, colNumber]
  data[,colNumber] <- NULL
  return (data)
} 
####################################################################
#             Load the unique information about data               #
####################################################################
getData <- function(pathwayData) {
# The rawData contains the probe names and all of the measurement 
# for each sample in the study. The pathway and the seperator for the
# data file is given. The data in this file will be saved as rawData.
  rawData <- read.delim(pathwayData, sep="\t")
  rawData <- rawData[order(colnames(rawData), decreasing = F )]
  rawData <- setRowname(rawData)
  return(rawData)
}

getTarget <- function(pathwatTargets) {
  # The information of the rawData is stored within a different file.
  # This file contains information like: the sample_nr, the condition,
  # the number of the plate and so on. The information of this file
  # is saved as targets.
  targets <- read.table(pathwatTargets, sep = ",", header = T)
  targets <- targets[order(targets$Sample),]
  return(targets)
}

####################################################################
#                         Filter functions                         #
####################################################################
filterGenes <- function(tabel, fdr, logfc) {
  # This function is made to calculate the unique genes for the 
  # differential expressions. The genes with a FDR below a number of 
  #0.05 are returned. To adjust the amount of returned genes, an 
  # adjusted FDR and logFC can be given. 
  # The if-else statements makes sure that the fdr and the logfc
  # values are optional.
  if (missing(fdr) & missing(logfc)) {
    uniqueGenes <- which(tabel[[1]]$FDR < 0.05)
  } else if (missing(logfc)) {
    uniqueGenes <- which(tabel[[1]]$FDR < fdr)
  } else if (missing(fdr)) {
    uniqueGenes <- which(abs(tabel[[1]]$logFC > logfc))
  } else {
    uniqueGenes <- which(abs(tabel[[1]]$logFC > logfc) & tabel[[1]]$FDR < fdr)
  }
  uniqueGenesList <- tabel[uniqueGenes,]
  BioM.uniqueGenes <- BioM[uniqueGenes,]
  data <- list(uniqueGenesList, BioM.uniqueGenes)
  return (data)
}

filterResult <- function(tabel, bioMList, calculation) {
# The function filterResults gets the genes and their 
# information according to the given calculation
# The BioM list is adjusted immediately as well.
  uniqueGenes <- calculation
  uniqueGenesList <- data.frame(tabel[uniqueGenes,])
  BioM.uniqueGenes <- bioMList[uniqueGenes, ]
  rownames(uniqueGenesList) <- BioM.uniqueGenes[,3]
  return(uniqueGenesList)
}

####################################################################
#                          Plot function                           #
####################################################################
mdsPlot <- function(pathway){
# The mdsPlot functions makes sure that a MDS plot is made
# The labels show the condition of the samples.
  pdf(pathway) 
  plotMDS.DGEList(dge, col=col_cell_age, main="MDS plot", labels = targets$Conditie)
  dev.off() 
}

plotMostExpr <- function(results, calc1, calc2, amount, path, type) {
# The function plotMostExpr makes sure that a heatmap is made
# of the genes that are returned by the function filterResult.
# There are two heatmaps in stead of two, these two conditions
# are found the most informative.
  results1 <- filterResult(results[[2]], results[[5]], calc1)
  results2 <- filterResult(results[[3]], results[[6]], calc2)
  results3 <- filterResult(results[[1]], results[[4]], resultsOlderMice[[1]][[1]]$FDR <  1E-8)
  write.csv(results1, paste(path, "ResultsForOneHeatmap", type, "MA.txt", sep = ""), row.names = rownames(results1), eol=",\n", quote = F)
  write.csv(results3, paste(path, "ResultsForOneHeatmap", type, "MG.txt", sep = ""), row.names = row.names(results3), eol=",\n", quote = F)
  write.csv(results2, paste(path, "ResultsForOneHeatmap", type, "IE.txt", sep = ""), row.names = row.names(results2), eol=",\n", quote = F)
  pdf(paste(path, "AgeAndLinear", type, ".pdf", sep=""))
  heatmap.2(M3[match(rownames(results1), rownames(M3)), amount], ColSideColors= col_cell_age[amount], col=hmcol, cexRow = 1, trace = "none", scale = "row", main="Main effect age")
  heatmap.2(M3[match(rownames(results2), rownames(M3)), amount], ColSideColors= col_cell_age[amount], col=hmcol, cexRow = 1, trace = "none", scale = "row", main="Interaction Effect")
  dev.off()
}

pcaPlot <- function(pathway) {
  ####################################################################
  #                       Preparing data for PCA                     #
  ####################################################################
  # The dataset is rotated and after this the coloms will be filtered
  # When the sum of the data is a negative number, this colomn will be deleted.
  tset <- t(M2)
  tset <- tset[,apply(tset, 2, sum) > 0]
  pca <- prcomp(tset, cor=T, scores= T, scale=T) 
  ####################################################################
  #                           Creating PCA                           #
  ####################################################################
  # open3d is used to save the data
  open3d()
  # The PCA plot is made with the use of the amount of samples.
  plot3d(pca$x, col=col_cell_age, size = "5")
  # rgl.postscript writes the current figure into a pdf (so you can rotate the plot before saving). 
  rgl.postscript(paste(pathway, "Plots/PCA.pdf", sep=""), "pdf")
}

makeVulcanoPlot <- function(table, result, name) {
  # Creates the vulano plot with the use of the results gained
  # from the function calFDRandLogFC.
  plot(as.data.frame(table[[1]][,1])[[1]], -log(as.data.frame(changeFDR(table[[1]][,5]))[[1]], 10), main=name, pch=20, xlab="log2FC)", ylab="-log10 (FDR p-value)")
  points(as.data.frame(result[[1]][,1])[[1]], -log(as.data.frame(changeFDR(result[[1]][,5]))[[1]], 10), cex = 1, col= "darkred", pch= 20)
  points(as.data.frame(result[[2]][,1])[[1]], -log(as.data.frame(changeFDR(result[[2]][,5]))[[1]], 10), cex = 1, col= "red", pch= 20)
  text(as.data.frame(result[[3]][,1])[[1]], -log(as.data.frame(changeFDR(result[[3]][,5]))[[1]], 10), result[[4]], cex = 0.75, col= "black")
}
####################################################################
#                        Calculate functions                       #
####################################################################
calFDRandLogFC <- function (tablename, pval1, pval2, pval3, logfc, logfc2) {
  # Three different toptables are made to create the vulcano plot.
  # This depends on the values given when calling this function.
  # The first value shows the values below a certain pval. The two
  # other tables need a pval and a logfc to make sure that the found genes 
  # are more significant then the one above.
  Toptable1 <- tablename[abs(as.data.frame(tablename[,5])) < pval1,] 
  Toptable2 <- tablename[abs(as.data.frame(tablename[,5])) < pval2 & abs(as.data.frame(tablename[,1])) > logfc,]
  Toptable3 <- tablename[abs(as.data.frame(tablename[,5])) < pval3 & abs(as.data.frame(tablename[,1])) > logfc2,]
  Genes1 <- BioM[abs(as.data.frame(tablename[,5])) < pval3 &  abs(as.data.frame(tablename[,1])) > logfc2, 3]
  results <- list(Toptable1, Toptable2, Toptable3, Genes1)
  return (results)
}

calculateMaineffectsInteraction <- function(dge, design, pathwayDoc, pathwayPlot, length) {
  # This function makes sure that the main effect of the ages, the main effect of the genotype and the interaction model
  # of each dataset is calculated.
  # The dge and the design are necessary items to calculate the these effects.
  # glmFit conducts statistical tests to fit the negative binomial generalized linear model.
  # The function glmLRT, tests the likelihood of the coefficients.
  # topTags orders the data by the ranking of the p value ord the logFC.
  # The lasts steps contains the filtering of the genes with a FDR < 0.05, the visualisation of these genes and 
  # saving these unique gene names into a file.
  fit <- glmFit(dge, design)
  # With the use of different coefs, different genewise statistical tests are made.
  # The first one compares the WT against APP (main effect of genotype)
  lrt1 <- glmLRT(fit, coef=2)
  # The second comparison checks the main effect of age
  lrt2 <- glmLRT(fit, coef=3)
  # The third and last comparison checks the linear interaction.
  lrt3 <- glmLRT(fit, coef=4)
  
  # The 3 different statistical tests are extracted in a dataframe
  # These results will be checked within the own-made function filtergenes to check if they meet the requirements.
  Toptable1 <- topTags(lrt1, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  Toptable2 <- topTags(lrt2, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  Toptable3 <- topTags(lrt3, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  # The differential expressed genes between the two genotypes
  Toptable1.results <- filterGenes(Toptable1)
  geneInfo1 <- Toptable1.results[[2]]
  Toptable1.results <- Toptable1.results[[1]]
  # The differential expressed genes over time
  Toptable2.results <- filterGenes(Toptable2)
  geneInfo2 <- Toptable2.results[[2]]
  Toptable2.results <- Toptable2.results[[1]]
  # The differential expressed genes in the interaction effect genotype*time
  Toptable3.results <- filterGenes(Toptable3)
  geneInfo3 <- Toptable3.results[[2]]
  Toptable3.results <- Toptable3.results[[1]]
  
  # The last step is to plot the genes with their information in a heatmap.
  pdf(pathwayPlot) 
  heatmap.2(M2[match(rownames(Toptable1.results), rownames(M2)), length], ColSideColors= col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Main effect genotype")
  heatmap.2(M2[match(rownames(Toptable2.results), rownames(M2)), length], ColSideColors = col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Main effect age")
  heatmap.2(M2[match(rownames(Toptable3.results), rownames(M2)), length], ColSideColors = col_cell_age[length], cexRow = 0.01, trace = "none", scale = "row", main="Interaction effect")
  dev.off()
  # The unique genes are saved within a table so these gene names can be used for futher analysis. 
  write.table(rownames(Toptable1.results), paste(pathwayDoc, "main_genotype_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)
  write.table(rownames(Toptable2.results), paste(pathwayDoc, "main_age_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)
  write.table(rownames(Toptable3.results), paste(pathwayDoc, "interaction_result.txt", sep = ""), row.names = F, col.names=F, eol=",\n", quote = F)
  # Creation of a DE file of the main and the interaction effects.
  DE_Expression <- cbind(rownames(Toptable1[[1]]), Toptable1[[1]]$logFC, Toptable1[[1]]$FDR, Toptable2[[1]]$logFC, Toptable2[[1]]$FDR, Toptable3[[1]]$logFC, Toptable3[[1]]$FDR, BioM[,3:4])
  write.table(DE_Expression , paste(pathwayDoc, "LinearTimeDE.txt", sep=""), row.names = F,  col.names = c("Genes", "Main-effect Genotype (logFC)", 
                                                                                                           "Main-effect Genotype (FDR)", "Main-effect Age (logFC)", "Main-effect Age (FDR)", "Interaction Effect (logFC)", 
                                                                                                           "Interaction Effect (FDR)", "Gene Symbol", "Gene Description"), sep = "\t")
  data <- list(Toptable1.results, Toptable2.results, Toptable3.results, geneInfo1,geneInfo2, geneInfo3)
}

####################################################################
#                          Other functions                         #
####################################################################
changeFDR <- function(table) {
  # If the pvalue is below the 1E-70, the value will be adjucted 
  # to a value of 1E-70, to make sure that the image looks more clear.
  pval <- table[[1]][1]
  pval[pval < 1E-70,] <- 1E-70 
  return (pval)
}

createToptableResults <- function(lrt, pathway, fdr, pval){
  # Is there is no pathway given when calling this function, 
  # a pathway is used that is defined in the function itself.
  if (missing(pathway)) {
    pathway = paste(resultPathway, "Made_Documents/DE/", sep="")
  }
  # Data is stored within the toptables for further anaylsis.
  toptable <- topTags(lrt, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  # The results from the toptables are checked with the use of the own-made filtergenes function
  # This function returns the data which contain a FDR value below the 0.05
  if (missing(fdr)) {
    toptable.results <- filterGenes(toptable, pval)
  } else if (missing(pval)) {
    toptable.results <- filterGenes(toptable, fdr)
  } else {
    toptable.results <- filterGenes(toptable, fdr, pval)
  }
  toptable.results <- toptable.results[[1]]
  # The comparison of the samples will be split at the ")" sign for the first time,
  # the second time the 2th string is split at the " -".
  # this results two genotype and age sets.
  # These information will be added together to create the filename.
  firstSplit <- strsplit(lrt$comparison, split = ")")
  secondSplit <- strsplit(firstSplit[[1]][2], split = " ")
  write.table(rownames(toptable.results), paste(pathway, "DE_", firstSplit[[1]][3], "_vs_", secondSplit[[1]][1], ".txt", sep = ""), eol=",\n", quote = F, row.names = F, col.names = F)
  output <- list(toptable, toptable.results)
  return (output)
} 

saveInfoDE <- function(result, fileName1, fileName2, fileName3) {
# The Differential expression genes are saved into the vector, the 
# second step is to create a list with information that can be used
# in the columns for the differential expression data. The data of 
# the DE data set and the list of names will be used to create a 
# text file with all of the information.
  DE.ExpressionOM_mainGenotype <- cbind(rownames(result[[1]][[1]]), result[[1]][[1]]$logFC, result[[1]][[1]]$FDR, result[[4]][,3:4])
  DE.ExpressionOM_mainAge <- cbind(rownames(result[[2]][[1]]), result[[2]][[1]]$logFC, result[[2]][[1]]$FDR, result[[5]][,3:4])
  DE.ExpressionOM_interaction <- cbind(rownames(result[[3]][[1]]), result[[3]][[1]]$logFC, result[[3]][[1]]$FDR, result[[6]][,3:4])
  geneColsOM_mainGenotype <- c("Genes", "logFC: Main-effect Genotype","FDR: Main-effect Genotype", "Gene Symbol", "Gene Description")
  geneColsOM_mainAge <- c("Genes", "logFC: Main-effect Age","FDR: Main-effect Age", "Gene Symbol", "Gene Description")
  geneColsOM_mainInteraction <- c("Genes", "logFC: Interaction Effect Age:Genotype",  "FDR:  Interaction Effect Age:Genotype", "Gene Symbol", "Gene Description")
  write.table(DE.ExpressionOM_mainGenotype, paste(resultPathway, "Made_Documents/DE_Files/", fileName1, sep = ""), row.names = F, col.names = geneColsOM_mainGenotype, sep = "\t")
  write.table(DE.ExpressionOM_mainAge, paste(resultPathway, "Made_Documents/DE_Files/", fileName2, sep = ""), row.names = F, col.names = geneColsOM_mainAge, sep = "\t")
  write.table(DE.ExpressionOM_interaction, paste(resultPathway, "Made_Documents/DE_Files/", fileName3, sep = ""), row.names = F, col.names = geneColsOM_mainInteraction, sep = "\t")
}