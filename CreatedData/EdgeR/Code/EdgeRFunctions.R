####################################################################
# Author    : M. Dubbelaar
# Date      : 10-sept-2015
# File Name : loadingAppFile.R
# Purpose   : Loads the data from the APP23 research and the information of the samples.
####################################################################
#                     Get all of the probe names                   #
####################################################################
setRowname <- function(data) {
  colNumber <- grep("probe", colnames(data))
  row.names(data) <- data[, colNumber]
  data[,colNumber] <- NULL
  return (data)
} 
####################################################################
#             Load the unique information about data               #
####################################################################
getData <- function(pathwayData, seperator) {
#  The rawData contains the probe names and all of the measurement for each sample in the study
  rawData <- read.delim(pathwayData, sep=seperator)
  rawData <- rawData[order(colnames(rawData), decreasing = F )]
  rawData <- setRowname(rawData)
  return(rawData)
}

getTarget <- function(pathwatTargets) {
  # The information of the rawData is stored within a different file.
  # This file contains information like: the sample_nr, the condition, the number of the plate and so on.
  targets <- read.table(pathwatTargets, sep = ",", header = T)
  targets <- targets[order(targets$Sample),]
  return(targets)
}

####################################################################
#                         Filter functions                         #
####################################################################
filterData <- function() {
  ####################################################################
  #                      Differential Expression                     #
  ####################################################################
  # A dge list is creates, this list is used to estimate the GLM 
  # common, trended and tagwise dispersion. This object holds the 
  # dataset that needs to be analysed by EgdeR and the calculations 
  # which where performed on the dataset
  dge <- DGEList(counts=rawData, group=factor(targets$Conditie) )
  ####################################################################
  #                      Loading Ensembl data                        #
  ####################################################################
  source("/home/mdubbelaar/APP23/CreatedData/loadingEnsemblData.R")
  ####################################################################
  #                       Estimation Dispersions                     #
  ####################################################################
  # The data is filtered with a cut-off of 1 count per million (cpm).
  # The data with 2 replicated samples will be saves into the vector.
  isExpr <- rowSums(cpm(dge)>1) >= 2
  # These expressed genes will be saved within the dge and BioM dataset.
  dge <- dge[isExpr, ]
  BioM <- BioM[isExpr,]
  
  # The raw library sizes are calculated with the use of calNormFactors
  # This is a step before calculating the estimates.
  dge <- calcNormFactors(dge)
  colnames(dge$counts) <- targets$Samples
  M2 <- cpm(dge, log=TRUE)
  M3 <- M2
  rownames(M3) <- BioM[,3]
}

filterGenes <- function(tabel, fdr, logfc) {
  # This function is made to calculate the unique genes for the differential expressions.
  # The genes with a FDR below a number of 0.05 are returned.
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
  pdf(pathway) 
  plotMDS.DGEList(dge, col=col_cell_age, main="MDS plot", labels = targets$Conditie)
  dev.off() 
}

plotMostExpr <- function(results, calc1, calc2, amount, path) {
  results1 <- filterResult(results[[2]], results[[5]], calc1)
  results2 <- filterResult(results[[3]], results[[6]], calc2)
  pdf(path)
  heatmap.2(M3[match(rownames(results1), rownames(M3)), amount], ColSideColors= col_cell_age[amount], cexRow = 1, trace = "none", scale = "row", main="Main effect age")
  heatmap.2(M3[match(rownames(results2), rownames(M3)), amount], ColSideColors= col_cell_age[amount], cexRow = 1, trace = "none", scale = "row", main="Linear Effect")
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
####################################################################
#                    Other (bigger functions                       #
####################################################################
createToptableResults <- function(lrt, pathway, fdr, pval){
  if (missing(pathway)) {
    pathway = "/home/mdubbelaar/Desktop/Results/APP23/Made_Documents/DE/"
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

saveInfoDE <- function(name, result, fileName1, fileName2, fileName3) {
  DE.ExpressionOM_mainGenotype <- cbind(rownames(result[[1]][[1]]), result[[1]][[1]]$logFC, result[[1]][[1]]$FDR, result[[4]][,3:4])
  DE.ExpressionOM_mainAge <- cbind(rownames(result[[2]][[1]]), result[[2]][[1]]$logFC, result[[2]][[1]]$FDR, result[[5]][,3:4])
  DE.ExpressionOM_Linear <- cbind(rownames(result[[3]][[1]]), result[[3]][[1]]$logFC, result[[3]][[1]]$FDR, result[[6]][,3:4])
  geneColsOM_mainGenotype <- c("Genes", "logFC: Main-effect Genotype","FDR: Main-effect Genotype", "Gene Symbol", "Gene Description")
  geneColsOM_mainAge <- c("Genes", "logFC: Main-effect Age","FDR: Main-effect Age", "Gene Symbol", "Gene Description")
  geneColsOM_mainLinear <- c("Genes", "logFC: Linear Effect Age:Genotype",  "FDR:  Linear Effect Age:Genotype", "Gene Symbol", "Gene Description")
  write.table(DE.ExpressionOM_mainGenotype ,paste("/home/mdubbelaar/Desktop/Results/", name,"/Made_Documents/DE_Files/", fileName1, sep = ""), row.names = F, col.names = geneColsOM_mainGenotype, sep = "\t")
  write.table(DE.ExpressionOM_mainAge , paste("/home/mdubbelaar/Desktop/Results/", name,"/Made_Documents/DE_Files/", fileName2, sep = ""), row.names = F, col.names = geneColsOM_mainAge, sep = "\t")
  write.table(DE.ExpressionOM_Linear , paste("/home/mdubbelaar/Desktop/Results/", name,"/Made_Documents/DE_Files/", fileName3, sep = ""), row.names = F, col.names = geneColsOM_mainLinear, sep = "\t")
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
  
  Toptable1.results <- filterGenes(Toptable1)
  geneInfo1 <- Toptable1.results[[2]]
  Toptable1.results <- Toptable1.results[[1]]
  
  Toptable2.results <- filterGenes(Toptable2)
  geneInfo2 <- Toptable2.results[[2]]
  Toptable2.results <- Toptable2.results[[1]]
  
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
  
  DE_Expression <- cbind(rownames(Toptable1[[1]]), Toptable1[[1]]$logFC, Toptable1[[1]]$FDR, Toptable2[[1]]$logFC, Toptable2[[1]]$FDR, Toptable3[[1]]$logFC, Toptable3[[1]]$FDR, BioM[,3:4])
  write.table(DE_Expression , paste(pathwayDoc, "LinearTimeDE.txt", sep=""), row.names = F,  col.names = c("Genes", "Main-effect Genotype (logFC)", 
                                                                                                           "Main-effect Genotype (FDR)", "Main-effect Age (logFC)", "Main-effect Age (FDR)", "Linear Effect (logFC)", 
                                                                                                           "Linear Effect (FDR)", "Gene Symbol", "Gene Description"), sep = "\t")
  data <- list(Toptable1.results, Toptable2.results, Toptable3.results, geneInfo1,geneInfo2, geneInfo3)
}

differentialExpression <- function(name, pathway) {   
  ####################################################################
  #                            MDS Plot                              #
  ####################################################################
  # The MDS plot is made and saved within a file with the use of pdf()
  # dev.off() is necessary to close the connection to the document.
  #mdsPlot("/home/mdubbelaar/Desktop/Results/Human/MDS.pdf")
  mdsPlot(paste(pathway, "Plots/MDS.pdf", sep=""))
  if (name == "APP23") {
    ####################################################################
    #                      Data Check (APP23 only)                     #
    ####################################################################
    CPMmatrix <- as.matrix(cpm(dge, log = T))
    rownames(CPMmatrix) <- BioM[,3]
    
    appGene <- which(rownames(CPMmatrix)=="App")
    appExpression <- CPMmatrix[appGene,]
    pdf(paste(pathway, "Plots/AppCheck.pdf", sep=""))
    #pdf("/home/mdubbelaar/Desktop/APP23_results/EdgeR/Plots/AppCheck.pdf")
    plot(appExpression, col=col_cell_age, pch=20, cex=3, ylab = "Expression of App")
    dev.off()
    
    source("/home/mdubbelaar/APP23/CreatedData/EdgeR/Code/EdgeRLinearTimeAPP23.R")
    source("/home/mdubbelaar/APP23/CreatedData/EdgeR/Code/EdgeRDispersionsAPP23.R")
    source("/home/mdubbelaar/APP23/CreatedData/saveDEData.R")
  } else if (name== "Human") {
    source("/home/mdubbelaar/APP23/CreatedData/EdgeR/Code/EdgeRDispersionsHuman.R")
  } else if (name == "CKp25") {
    source("/home/mdubbelaar/APP23/CreatedData/EdgeR/Code/EdgeRDispersionsCKp25.R")
  } else {
    print("Dataset is not known")
  }
  ####################################################################
  #                          Quality plots                           #
  ####################################################################
  # Plots the biological coefficient of variantion agains the log2 
  # counts per million. The common, trended ans tagwise BCV estimates
  # are shown in the plot
  plotBCV(dge, cex=0.4)
  # Plots the sample relations. The distances of the RNA-seq is 
  # calulated, the distances represent the coefficient of variation
  # of expression among the sampels.
  plotMDS.DGEList(dge, col=col_cell_age)
  ####################################################################
  #                    Spearman correlatie plot                      #
  ####################################################################
  # A heatmap with the correlation of the different genes are shown.
  # This correlation is calculated from the M2 dataset.
  cor <- cor(M2, method="spearman")
  rownames(cor) <- targets$Samples
  colnames(cor) <- targets$Conditie
  #pdf("/home/mdubbelaar/Desktop/Human_results/Plots/Spearman_Cor_Plot.pdf")
  pdf(paste(pathway, "Plots/Spearman_Cor_Plot.pdf", sep=""))
  heatmap.2(cor, symm = TRUE, col=greenred(350),
            trace="none", cexRow = 1 , cexCol = 1, ColSideColors= col_cell_age, RowSideColors=col_cell_age)
  dev.off()
  ####################################################################
  #                           PCA Plot                               #
  ####################################################################
  pcaPlot(pathway)
}