source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/Desktop/Onderzoek-APP23_RNASEQ/Code_APP23_RNASEQ/")
#biocLite("limma" )
#biocLite("biomaRt")
#biocLite("gplots")
#biocLite("edgeR")
library(biomaRt)
library(gplots)
library(limma)
library(edgeR)
library(statmod)
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("loading_file.r")
M1 <- APP23_data
row.names(M1) <- APP23_data[,1]
source("loading_ensembl_data.r")
dge <- DGEList(counts=M1[,2:25])

#A <- rowSums(dge$counts)
#isexp <- A > 50
#dge <- dge[isexp, keep.lib.size=F]

isExpr <- rowSums(cpm(dge)>1) >= 2
dge <- dge[isExpr, keep.lib.sizes=F]
BioM <- BioM[isExpr,]
dge <- calcNormFactors(dge)
rownames(dge) <- BioM[,3]

# scale normalization based on the wildtypes and the APP-mice
object <- substring(targets$Conditie,1,1)
plotMDS(dge, labels = object, top = 50, col=ifelse(object=="W", "blue", "darkgoldenrod"))

design <- model.matrix(~factor(targets$Conditie), data = dge$samples)
####################################################################
#                      Loading Ensembl data                        #
####################################################################
#source("loading_ensembl_data.r")

# The voom transformation, produces an EList object, by using the design matrix.
v <- voom(dge, design, plot = T)

fit <- lmFit(v, design)
fit <- eBayes(fit)
topTable(fit, coef = ncol(design))

plotMDS(dge, labels = targets$Conditie, col=as.numeric(targets$Conditie))
legend("topright", legend = unique(targets$Conditie), col=unique(targets$Conditie), pch=15, cex=.6)

####################################################################
#                            Visualization                         #
####################################################################
# Hieronder kan een matrix gemaakt worden (probeer deze hetzelfde te maken als bij de manier van EdgeR)
cont.matrix <- makeContrasts( 
  All.minus6_8weeks = (-WT.24-WT.18-WT.12)/3-(-HET.24-HET.18-HET.12)/3,
  WT.HET12_18 = (-WT.18-WT.12)/2-(-HET.18-HET.12)/2,
  WT.HET18_24 = (-WT.24-WT.18)/2-(-HET.24-HET.18)/2,
  WT18_24 = WT.24-WT.18,
  HET18_24 = HET.24-HET.18,
  All.minus12months = (-WT.24-WT.18-WT.06_08)/3-(-HET.24-HET.18-HET.06_08)/3,
  WT.HET6.8_18 = (-WT.18-WT.06_08)/2-(-HET.18-HET.06_08)/2,
  WT.HET6.8_24 = (-WT.24-WT.06_08)/2-(-HET.24-HET.06_08)/2,
  All =(-WT.24-WT.18-WT.12-WT.06_08)/4-(-HET.24-HET.18-HET.12-HET.06_08)/4,
  levels=targets$Conditie
  )
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
results <- decideTests(fit2, method = "global")

#vennDiagram(results[,1:3], cex = 0.4, circle.col = c(1:3))
vennDiagram(results[,c(1,2,3,9)], cex = 0.4)
laterStageGene <- which(results[,9] & results[,3] & results[,2] & results[,1])
laterStageGene <- results[laterStageGene,]
rownames(laterStageGene)

vennDiagram(results[,c(6,7,8,9)], cex = 0.4)
earlierStageGene <- which(results[,9] & results[,8] & results[,6] & results[,7])
earlierStageGene <- results[earlierStageGene,]
rownames(earlierStageGene)

vennDiagram(results[,c(4,5,9)], cex = 0.4)
wt.against.het <- which(results[,5] & results[,4] & results[,9])
wt.against.het <- results[wt.against.het,]
rownames(wt.against.het)

summary(results)
heatmap.2(cont.matrix, cexRow = 0.8, cexCol = 0.5 )