source("http://bioconductor.org/biocLite.R") 
setwd("/home/mdubbelaar/APP23/RScripts/WGCNA/")
####################################################################
#install.packages("WGCNA")
#biocLite("GO.db")
#biocLite("preprocessCore")
#biocLite("biomaRt")
#biocLite("flashClust")
#biocLite("org.Mm.eg.db")
#biocLite("WGCNA")
#biocLite("edgeR")
#biocLite("gplots")
#install.packages("corrplot")
####################################################################
library("corrplot")
library("WGCNA")
library("cluster")
library("biomaRt")
library("edgeR")
library("flashClust")
library("org.Mm.eg.db")
library("gplots")
allowWGCNAThreads()
ALLOW_WGCNA_THREADS=8
resultPathway <- "/home/mdubbelaar/Desktop/APP23_results/WGCNA/"
####################################################################
#                            Loading data                          #
####################################################################
source("../EdgeR/Code//EdgeRFunctions.R")
M1 <- getData("/media/mdubbelaar/6CEC0BDEEC0BA186/1507_Holtman_RNAseq/run01/results/expression/expressionTable/expression_table02.genelevel.GRCm38.v76.htseq.txt.table")
targets <- getTarget("../../Targets.csv")
source("../loadingEnsemblData.R")
hmcol <- colorRampPalette(c("purple", "lightgrey", "green")) (n=99)
####################################################################
#                          Filtering data                          #
####################################################################
# A DGEList is made, this is an object that contains the counts, a
# group indicator for each column, the library size and the feature
# annotation.
dge <- DGEList(M1, group=factor(targets$Conditie))

CPMmatrix <- as.matrix(cpm(dge, log = T))
BioM_CPMmatrix <- BioM[match(rownames(CPMmatrix), BioM[,1]),]
rownames(CPMmatrix) <- BioM_CPMmatrix[,3]
appGene <- which(rownames(CPMmatrix)=="App")
appExpression <- CPMmatrix[appGene,]

# The function calcnormfactors makes sure that the raw library sizes
# are scaled by calculation the normalization factors.
dge <- calcNormFactors(dge, "upperquartile")
# The design for the data is given, the condition shows the genotype
# and the age of each sample.
design <- model.matrix(~0+factor(targets$Conditie), data = dge$samples)  

# The genes will be filtered
# The counts of each gene within each sample will be divided by the library
# size of each sample. To make sure that the number is not to low the result 
# will be multiplied by 1e6. 
m <- 1e6 * t(t(dge$counts) / dge$samples$lib.size)
# The sum of the row will be taken, the expression of each gene is 
# followed in this step over all samples. All genes that have a value
# higher than 2 are saved and used further.
goodGenes <- rowSums(m > 1) >= 2
dge <- dge[goodGenes,]
BioM <- BioM[goodGenes,]

# The following steps are used to estimate the dispersions within the
# dataset, this estimates will be used further.
D <- estimateGLMCommonDisp(dge, design)
D <- estimateGLMTrendedDisp(D, design)
D <- estimateGLMTagwiseDisp(D, design)
# Generation of a normalized data set, the results are returned in the
# form of the log2 values.
M2 <- cpm(D, log=TRUE)

####################################################################
#                      Soft Threshold beta                         #
####################################################################
# The code below makes sure that the data is filtered
# The expression of each gene above the 25% Quantile is saved
# and used again in the dataset.
INCL <- apply(M2, 1, sd) > quantile(apply(M2, 1, sd))[2]
M2 <- M2[INCL,]
BioM <- BioM[INCL,]

powers <- c(1:20)
# The function pickSoftThreshold is used to analyse the scale free
# topology with the use of multiple soft-thresholding powers.
# The goal is to help the user pick an apporpriate power for the network 
# construction. The beta value 5 will be used, the R^2 that is reached will be 0.862
sft2 <- pickSoftThreshold(t(M2),powerVector=powers)
plot(sft2$fitIndices[,1],-sign(sft2$fitIndices[,3])*sft2$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="SFT2, signed R^2",type="n",main=paste("Scale independence"))
text(sft2$fitIndices[,1],-sign(sft2$fitIndices[,3])*sft2$fitIndices[,2],
     labels=powers,col="black")
abline(h=0.85, col="red")
plot(sft2$fitIndices[,1],sft2$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft2$fitIndices[,1],sft2$fitIndices[,5],labels=powers,col="red")
###################################################################
#                     Stepwise module detection                   #
###################################################################
# You need to double the power when the soft power is used.
# The reason for this lies in the signed and unsigned type.
# The unsigned type gets a softpower 5 (which is the power that is
# obtained before), the signed doubles the amount of power.
softPower <- 10
adjacencyM2 <- adjacency(t(M2),power=softPower,type="signed");
# The diagonal value is put to 0, this will make sure that the 
# gene doens't have a relationship with itself.
diag(adjacencyM2) <- 0

# This function and the function that is commented are the same
# These function make sure that the distance among the different 
# genes are calculated.
dissTOMM2   <- 1-TOMsimilarity(adjacencyM2, TOMType="signed")
# dissTOMM2   <- TOMdist(adjacencyM2, TOMType="signed")
# The function flashClust performs a hierarchical clustering on the data.
geneTreeM2  <- flashClust(as.dist(dissTOMM2), method="average")

# cutreeDynamic makes sure that the data from dissTOMM2 and geneTreeM2
# can be visualized together.
dynamicMods <- cutreeDynamic(geneTreeM2, distM = dissTOMM2, pamStage= T, 
                             deepSplit = 1, pamRespectsDendro = T,
                             minClusterSize = 100, cutHeight = 0.9999)
# The different clusters that are found in the dynamicModel can be 
# used by labels2colors to get the colors of the different modules that
# are observed and to use these models further.
dynamicCols <- labels2colors(dynamicMods)
# The creation of the dendrogram and the observed color modules can be
# seen in the plot that is generated by plotDendroAndColors.
plotDendroAndColors(geneTreeM2, dynamicCols, "Dynamic Tree Cut", dendroLabels = F,
                    hang = 0.03, addGuide = T, guideHang = 0.05, main="Gene dendrogram and module colors")

# The eigengenes are calculated, this eigengenes are used to calculate 
# the relation among the co-expression modules.
MEList <- moduleEigengenes(t(M2), colors = labels2colors(dynamicCols))
# The values of the eigengenes are stored into MEs
MEs <- MEList$eigengenes
# Calculating the dissimilarities among the module eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")
# The clustering of the different modules can be seen in the 
# plot that is made using the METree.
plot(METree, main="CLustering of eigengenes modules")
abline(h = "0.25", col="red")
# Colors that are too close in the gene expression network will be
# merged together. This will give a beter overview of the data.
merge <- mergeCloseModules(t(M2), dynamicCols, cutHeight = 0.25)
moduleCols <- merge$colors
mergedMEs <- merge$newMEs
# The new merged dendrogram and module colors can be seen in the 
# following plotDendroAndColors.
pdf(paste(resultPathway,"MergedDynamicTreeCut.pdf", sep=""))
plotDendroAndColors(geneTreeM2, moduleCols, "Dynamic Tree Cut",
                    dendroLabels = F, hang = 0.03, addGuide = T,
                    guideHang = 0.05)
dev.off()
###################################################################
#                    Relating Modules to traits                   #
###################################################################
# The Module eigengene are recalculated with the color labels
MEs0 <- moduleEigengenes(t(M2), moduleCols)$eigengenes
MEs <- orderMEs(MEs0)

# The Correlation and the pValue of the module eigengenes are calculated
# for some of the traits (Genotype, Age, development, APP aging and WT aging)
modTraitCorNormal <- cor(MEs, as.data.frame(targets[c(12, 7, 13:15)]), use="p")
modTraitPNormal <- corPvalueStudent(modTraitCorNormal, nrow(t(M2)))
# A text matrix is made, this matrix contains all of the 
# correlation values and the pvalues for each module and trait.
textMatrixNormal <- paste(signif(modTraitCorNormal,2), " \n(",
                          signif(modTraitPNormal, 1), ")", sep="")
dim(textMatrixNormal) <- dim(modTraitCorNormal[,1])

# The correlation values and the p-values are shown in the heatmap correlation plot.
pdf(paste(resultPathway,"Plots/Module-trait_RelationshipNormal.pdf", sep=""))
par(mar = c(6, 8.5, 3, 3))
# The rows of the heatmap correspond to the module eigengene and each row
# is known to a traits. The cells in the heatmap contain the correlation
# adn the p-value.
labeledHeatmap(Matrix = modTraitCorNormal, xLabels = c("Genotype","Age", "Development", "APP23 Aging", "WT Aging"), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = F, colors = hmcol, 
               textMatrix = textMatrixNormal, setStdMargins = F, cex.text = 0.8, zlim=c(-1, 1), cex.lab.x = 1.5,
               cex.lab.y = 0.0000001, main="Module-trait relationships")
dev.off()

# The Correlation and the pValue of the module eigengenes are calculated
# for some of the traits (AB, AB40 and AB42)
# The code below is without the log values among the ab values. 
modTraitCorAB <- cor(MEs, as.data.frame(targets[c(9:11)]), use="p")
modTraitPAB <- corPvalueStudent(modTraitCorAB, nrow(t(M2)))
# A text matrix is made, this matrix contains all of the 
# correlation values and the pvalues for each module and trait.
textMatrixAB <- paste(signif(modTraitCorAB,2), " \n(",
                      signif(modTraitPAB, 1), ")", sep="")
dim(textMatrixAB) <- dim(modTraitCorAB[,1])

# The correlation values and the p-values are shown in the heatmap correlation plot.
pdf(paste(resultPathway,"Module-trait_RelationshipAB.pdf", sep=""))
par(mar = c(6, 8.5, 3, 3))
# The rows of the heatmap correspond to the module eigengene and each row
# is known to a traits. The cells in the heatmap contain the correlation
# adn the p-value.
labeledHeatmap(Matrix = modTraitCorAB, xLabels = c("AB", "AB40", "AB42"), yLabels = names(MEs),
               ySymbols = names(MEs), colorLabels = F, colors = hmcol, 
               textMatrix = textMatrixAB, setStdMargins = F, cex.text = 0.8, zlim=c(-1, 1), cex.lab.x = 1,
               main="Module-trait relationships")
dev.off()
###################################################################
# The AppAndAB plot shows the different expressions of the APP gene
# expression and the AB, AB40 and AB42 ELISA values.
pdf(paste(resultPathway,"AppAndAB.pdf", sep=""))
source("..//plotColors.R")
par(mfrow=c(4,1), mar=c(1,2,1,0))
plot(appExpression, col=col_cell_age, pch=20, cex=3, main="APP")
plot(targets$Amyloid..g.wet.weight, col=col_cell_age, pch=20, cex=3, main="AB")
plot(targets$Amyloid40..g.wet.weight, col=col_cell_age, pch=20, cex=3, main="AB40")
plot(targets$Amyloid42..g.wet.weight, col=col_cell_age, pch=20, cex=3, main="AB42")
dev.off()

# The following plots show the correlation of the APP, AB, AB40 and
# AB42 all together. 
pdf(paste(resultPathway,"CorrelationAB.pdf", sep=""))
appData <- data.frame(appExpression, targets[9:11])
colnames(appData) <- c("APP", "AB", "AB40", "AB42")
corrplot(cor(appData), method = "number", diag = F)
dev.off()
###################################################################
# Get all of the modules names (also known as the colors)
modName <- substring(names(MEs), 3)
# The association of all genes are quantified with the traits by definint the 
# Gene siginificance (GS). This describes the absolute value of the correlation
# between the gene and trait. The quantitative measure of the module membership
# is also defines. These values can be used to quantify the similarity among the genes
geneModuleMembership <- as.data.frame(cor(t(M2), MEs, use = "p"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(t(M2))))
names(geneModuleMembership) <- paste("MM", modName, sep="")
names(MMPvalue) <- paste("p.MM", modName, sep = "")
geneTraitSignificance <- as.data.frame(cor(t(M2), targets[7], use = "p"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(t(M2))))
names(geneTraitSignificance) <- paste("GS.", names(targets[7]), sep = "")
names(GSPvalue) <- paste("p.GS.", names(targets[7]), sep = "")

# The GS and the MM values are used to identify genes with a high
# significance for the ab. 
# The code below here shows the scatterplot of all the modules.
pdf(paste(resultPathway,"ModuleMemberships.pdf", sep=""))
for (i in unique(modName)) {
  if (i != "grey") {
  module <- i
  column <- match(module, modName)
  moduleGenes <- moduleCols==module
  
  par(mfrow = c(1,1))
  verboseScatterplot(geneModuleMembership[moduleGenes, column],
                     geneTraitSignificance[moduleGenes, 1],
                     xlab = paste("Module MEmbership in", module, "module"),
                     ylab = "Gene significance for Age", main = paste(
                       "Module membership vs. gene significance\n"),
                     col = module)
  }
}
dev.off()
###################################################################
# The modules that contain a high association with the ab trait
# are identified with the use of the module membership measures.
# These results will be merges and will be saved into a csv file.
# The gene information that corresponds to the rownames in the 
# M2 data set are used.
BioM <- BioM[match(rownames(M2), BioM[,1]),]
# The geneInfo0 data frame contains gene information and the colors
# of the modules, together with the significance and the GS P-values.
geneInfo0 <- data.frame(ensemble=colnames(t(M2)), geneSymbol=BioM[,3],
                       moduleColor = moduleCols, geneTraitSignificance,
                       GSPvalue)
# The mod order can only take one column, so the column ab is chosen
# to be used in the forloop.
# The variables of the AB, the age and the genotype are used as a trait.
traits <- as.data.frame(targets[c(12, 7, 9:11)])
names(traits) <- c("Genotype", "Age","AB", "AB40", "AB42")
modOrder <- order(-abs(cor(MEs, traits[,2], use = "p")))
# The forloop check every cloumn in geneModuleMembership.
# Getting the information like the MM data and the p value of the 
# MM data of the corresponding module
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames <- names(geneInfo0)
  geneInfo0 <- data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
  MMPvalue[, modOrder[mod]])
  names(geneInfo0) <- c(oldNames, paste("MM.", modName[modOrder[mod]], sep = ""),
                        paste("p.MM.", modName[modOrder[mod]], sep = ""))
  }
# Orders the dataset
geneOrder <- order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Age))
geneInfo <- geneInfo0[geneOrder, ]
# This data will be written into the csv file.
write.csv(geneInfo, paste(resultPathway, "geneInfoAge.csv", sep=""))
###################################################################
#                    Visualization on networks                    #
###################################################################
# The MEs are ordered, this is to make sure that clusters can be
# seen in the correlation heatmap. The traits are also added to
# see the correlation of the traits among the different clusters.
MET <- orderMEs(cbind(MEs0, traits))
# A dendrogram and a heatmap is created to visualize the data. 
# It shows the clustering of the different modules and the 
# traits.
pdf(paste(resultPathway,"VisualizationNetwork.pdf", sep=""))
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), 
                      marHeatmap = c(3,4,1,2), cex.lab=0.8,
                      xLabelsAngle = 90)
dev.off()
###################################################################
#                 Relationship among the Modules                  #
###################################################################
pdf(paste(resultPathway,"ModulesBarplot.pdf", sep=""))
# This part of the code creates barplots that describe the expression of 
# each sample.
for (which.module in names(table(modName))){
  par(mfrow=c(1,1))
  # The column that can be used is determined by ME
  # since the ME is stored as ME[color]
  # This modules is used for the barplot.
  ME = MEs[, paste("ME",which.module, sep="")]
  barplot(ME, col=which.module, cex.main=2, 
          ylab="eigengene expression",xlab="array sample", main=which.module)}
dev.off()

pdf(paste(resultPathway,"ModulesBoxplot.pdf", sep=""))
# This part of the code create boxplots that describe the expression of 
# each condition.
for (i in 1:length(colnames(MEs))){
  MEs <- MEs[order(colnames(MEs))]
  # The names of the columns are identified and saved as name
  name <- colnames(MEs)[i]
  # The color of the modules is saved as color.
  color  <- unique(modName)[order(unique(modName))][i]
  # The boxplot with the expression of each condition is shown.
  # unfortunately the levels needed to be given by hand, to make sure that
  # it is not sorted alphabetically.
  verboseBoxplot(as.numeric(MEs[,name]), ordered(interaction(targets$Conditie), 
                 levels = c("WT.02", "HET.02", "WT.06", "HET.06", "WT.18", "HET.18", "WT.24", "HET.24")), 
                 main=name, las=2, xlab="", ylab="", notch=FALSE, col= color)
  }
dev.off()

# Is used to determine the amount of genes that are found in the given module.
#length(rownames(t(scale(t(M2))[, moduleCols == "red"])))

pdf(paste(resultPathway,"Heatmaps.pdf", sep=""))
for (which.module in names(table(modName))){
  # This part of the script shows the barplots again, but only to make sure
  # that the visualisation of the heatmap is easier to see.
  # The heatmap shows the underexpression (green) and overexpression (red)
  # of a given module for each sample.
  ME <- MEs[, paste("ME", which.module, sep = "")]
  par(mfrow=c(2,1), mar=c(1,4,1,1))
  barplot(ME, col=which.module, cex.main=1)
  par(mar=c(1,4.8,0.5,1.8))
  plotMat(t(scale(t(M2))[, moduleCols == which.module]), nrgcols = 30,
             rlabels = F)
}
dev.off()
###################################################################
#                         UserListEnrichment                      #
###################################################################
# The genModule from the geneInfoAge.csv is used to get the module
# of a specific gene.
geneModule <- read.csv(paste(resultPathway, "Modulemembership.csv", sep=""))
# The obtained data from geneModule is ordered by the ensemble id and
# saved into geneModuleData with the ensemble id, the gene name and the 
# module where the gene is categorized in.
geneModuleData <- geneModule[order(geneModule[,2]),2:4]
# The following step makes sure that the data is returned without NA's
filteredGeneModules <- na.omit(geneModuleData)
# Listgenes contains the unique genes in the filteredGeneModules dataset.
listGenes <- unique(as.character(filteredGeneModules[,2]))
# Categories contain all of the different modules.
categories <- filteredGeneModules[,3]
# The userListEnrichtment function is used to find interesting modules.
# For this project we used the available lists to compaire the current data.
# The results can be seen in the Enrichment.csv.
userListEnrichment(toupper(listGenes), labelR = categories,
                                  nameOut = paste(resultPathway, "userListEnrichment/Enrichment.csv", sep=""), 
                                  useBrainList=T, useBloodAtlases = T,
                                  useStemCellLists = T, useBrainRegionMarkers = T,
                                  useImmunePathwayLists = T, usePalazzoloWang = T, omitCategories = "grey")
###################################################################
#                  Network Visualization Software                 #
###################################################################
# The hub gene within each module can be foound with the use of 
# the function chooseTopHubInEachModule.
hubs <- chooseTopHubInEachModule(t(M2), moduleCols, omitColors = "grey",
                                 power=5, type = "signed")
# The ensemble id will be converted into the gene name after that.
BioM_hubs <- BioM[match(as.matrix(hubs)[,1], BioM[,1]),]
hubs <- cbind(hubs, BioM_hubs[,3])

# The idea of this part in the code it to make a file that can be
# used in VisANT. Each of the genes that are avaiable in the 
# M2 dataset are saved as probes.
probes <- colnames(t(M2))

# Each module is used in this function.
for (module in names(table(modName))){
  # The right module is selected.
  inModule <- moduleCols==module
  # The probes that are available for the module is saved as modProbes.
  modProbes <- probes[inModule]
  # The distances among the genes within a module are saved into modTOM
  modTOM <- dissTOMM2[inModule, inModule]
  dimnames(modTOM) <- list(modProbes, modProbes)
  # The function exportNetworkToVisant makes sure that the available data 
  # is saved into a txt file in a way that it can be used into VisANT.
  exportNetworkToVisANT(modTOM, file= paste(resultPathway, "VisOutput/VisANTInput-", module, ".txt", sep=""),
                        weighted = T, threshold = 0, probeToGene = data.frame(rownames(M2), BioM[,3]))
  
  # softConnectifity calculated the connectivity of each node that are found 
  # within one module. This data will be saved in to the dataframe IMConn
  IMConn <- softConnectivity(t(M2)[, modProbes])
  # The top 30 genes within this dataset are saved.
  top <- rank(-IMConn) <= 50
  # And a file that is saved.
  exportNetworkToVisANT(modTOM[top, top], file= paste(resultPathway, "VisOutput/VisANTInput-", module, "-top50.txt", sep=""),
                        weighted = T, threshold = 0, probeToGene = data.frame(rownames(M2), BioM[,3]))
}
# Directory is set to save the environment since some steps like
# calculation the TOM takes a lot of time.
save.image()
