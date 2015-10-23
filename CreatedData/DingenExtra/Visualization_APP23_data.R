####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : Visualization_APP23_data.R
# Purpose   : 
# Used Files: loading_file.r
#             loading_ensembl_data.r
#             estimating_dispersions.r
#             meer_testen.r
#             Creation_MDS_plots.r
#             Creating_PCA.r
####################################################################
#              Installing all of the necessary packages            #
####################################################################
source("http://bioconductor.org/biocLite.R")
setwd("/home/mdubbelaar/Desktop/Onderzoek-APP23_RNASEQ/Code_APP23_RNASEQ/")
library(limma)
library(edgeR)
library(biomaRt)
library(gplots)
####################################################################
#               Reading of the data and the targets                #
####################################################################
source("loading_file.r")
M1 <- APP23_data
row.names(M1) <- APP23_data[,1]
dge <- DGEList(counts=M1[,2:25], group=factor(targets$Conditie) )
####################################################################
#                      Loading Ensembl data                        #
####################################################################
source("loading_ensembl_data.r")
####################################################################
#                       Estimation Dispersions                     #
####################################################################
# The data is filtered with a cutoff of 1 count per million (cpm).
# And there must be at least two replicate samples in each group.
isExpr <- rowSums(cpm(dge)>1) >= 2
dge <- dge[isExpr, ]
BioM <- BioM[isExpr,]

dge <- calcNormFactors(dge)
design <- model.matrix(~0+factor(targets$Conditie), data = dge$samples)

colnames(design) <- levels(dge$samples$group)
M2 <- cpm(dge, log=TRUE)

source("estimating_dispersions.r")
source("meer_testen.r")
####################################################################
#                    Differential Expression                       #
####################################################################
# De DE_expression moet worden aangepast naar alle beschikbare toptables.
DE_expression <- cbind(rownames(Toptable1[[1]]), Toptable1[[1]]$logFC, Toptable1[[1]]$FDR, Toptable2[[1]]$logFC, Toptable2[[1]]$FDR, Toptable3[[1]]$logFC, Toptable3[[1]]$FDR, Toptable4[[1]]$logFC, Toptable4[[1]]$FDR, 
                       Toptable5[[1]]$logFC, Toptable5[[1]]$FDR, Toptable6[[1]]$logFC, Toptable6[[1]]$FDR, Toptable7[[1]]$logFC, Toptable7[[1]]$FDR, Toptable8[[1]]$logFC, Toptable8[[1]]$FDR, 
                       Toptable9[[1]]$logFC, Toptable9[[1]]$FDR, 
                       #Toptable10[[1]]$logFC, Toptable10[[1]]$FDR, Toptable11[[1]]$logFC, Toptable11[[1]]$FDR, Toptable12[[1]]$logFC, Toptable12[[1]]$FDR, 
                       #Toptable13[[1]]$logFC, Toptable13[[1]]$FDR, Toptable14[[1]]$logFC, Toptable14[[1]]$FDR, Toptable15[[1]]$logFC, Toptable15[[1]]$FDR, Toptable16[[1]]$logFC, Toptable16[[1]]$FDR, 
                       BioM[,2:4])
# Hieronder moeten de vergelijkingen worden genoteerd
colnames(DE_expression) <-c("Genes", "WT_VS_HET_logFC", "WT_VS_HET_FDR", "WT_12_18_24_VS_HET_12_18_24_logFC", "WT_12_18_24_VS_HET_12_18_24_FDR", "WT_12_18_VS_HET_12_18_logFC", "WT_12_18_VS_HET_12_18_FDR", "WT_18_24_VS_HET_18_24_logFC", "WT_18_24_VS_HET_18_24_FDR", "WT_18_VS_WT_24_logFC", "WT_18_VS_WT_24_FDR", "HET_18_VS_HET_24_logFC","HET_18_VS_HET_24_FDR", "WT_6/8_18_24_VS_HET_6/8_18_24_logFC", "WT_6/8_18_24_VS_HET_6/8_18_24_FDR", "WT_6/8_18_VS_HET_6/8_18_logFC", "WT_6/8_18_VS_HET_6/8_18_FDR", "WT_6/8_24_VS_HET_6/8_24_logFC", "WT_6/8_24_VS_HET_6/8_24_FDR", "ENSEMBL_Transcript","Gene_symbol", "Wiki-description")
#<- c("Genes", "HET-06/08_VS_WT-06/08_LogFC", "HET-06/08_VS_WT-06/08_FDR","HET-12_VS_WT-12_LogFC", "HET-12_VS_WT-12_FDR", "HET-18_VS_WT-18_LogFC", "HET-18_VS_WT-18_FDR", "HET-24_VS_HET-24_logFC", "HET-24_VS_WT-24_FDR",
#    "HET-06/08_VS_HET-12_LogFC", "HET-06/08_VS_HET-12_FDR", "HET-06/08_VS_HET-18_logFC", "HET-06/08_VS_HET-18_FDR", "HET-06/08_VS_HET-24_logFC", "HET-06/08_VS_HET-24_FDR", "HET-12_VS_HET-18_logFC", "HET-12_VS_HET-18_FDR", 
#    "HET-12_VS_HET-24_logFC", "HET-12_VS_HET-24_FDR",  "HET-18_VS_HET-24_logFC", "HET-18_VS_HET-24_FDR", "WT-06/08_VS_WT-12_logFC", "WT-06/08_VS_WT-12_FDR", "WT-06/08_VS_WT-18_logFC", "WT-06/08_VS_WT-18_FDR", 
#    "WT-06/08_VS_WT-24_logFC", "WT-06/08_VS_WT-24_FDR", "WT-12_VS_WT-18_logFC", "WT-12_VS_WT-18_FDR", "WT-12_VS_WT-24_logFC", "WT-12_VS_WT-24_FDR", "WT-18_VS_WT-24_logFC", "WT-18_VS_WT-24_FDR","ENSEMBL_Transcript","Gene_symbol", "Wiki-description")  

source("Creation_MDS_plots.r")
write.table(DE_expression, "../Plots_APP23_RNASEQ/Differential_expression.csv", row.names=FALSE, sep=":")
# Hieronder wordt gekeken naar de eerste 100 rows binnen de dataset van M2.
heatmap.2(M2[match(rownames(Toptable1[1:75,]), rownames(M2)),], ColSideColors = col_cell_age, density.info="none",  trace="none", scale = c("row"))
####################################################################
#                          Quality plots                           #
####################################################################
plotMeanVar(dge, show.raw.vars = T, NBline = T)
plotBCV(dge, cex=0.4)
plotMDS.DGEList(dge, col= col_cell_age)
####################################################################
#                    Spearman correlatie plot                      #
####################################################################
# Hieronder wordt er een heatmap aangemaakt die de samenhang van verschillende genen laat
# Deze samenhang is berekend na aanleiding van de correlatie uit de M2 dataset
cor <- cor(M2, method="spearman")
heatmap.2(cor, symm = TRUE, col=greenred(350),
          trace="none", cexRow = 1 , cexCol = 1, ColSideColors= col_cell_age, RowSideColors= col_cell_age)
####################################################################
#                           PCA Plot                               #
####################################################################
M2 <- cpm(dge, log=TRUE)
source("Creating_PCA.r")
####################################################################
#                            Heatmaps                              #
####################################################################
# bekijken wat hier gebeurd
# Hieronder wordt gekeken of er aan verschillende absolute waarden wordt gedaan
# Deze functie retourneerd TRUE of FALSE.
# De TRUE waarden worden opgeslagen in een tabel.
# Als laatste worden de namen aangepast naar de omschrijving van het gen.

WT6.8_18_24.vs.HET6.8_18_24.data <- abs(DE_expression[,14] < 1.5) & abs(DE_expression[,16] < 1.5) & abs(DE_expression[,18] < 1.5) & DE_expression[,15] <0.05 & DE_expression[,17] <0.05 & DE_expression[,19] <0.05
table(WT6.8_18_24.vs.HET6.8_18_24.data)
M2_WT6.8_18_24.vs.HET6.8_18_24 <- M2[WT6.8_18_24.vs.HET6.8_18_24.data,]
rownames(M2_WT6.8_18_24.vs.HET6.8_18_24) <- BioM[WT6.8_18_24.vs.HET6.8_18_24.data,3]

WT18_24.vs.HET18_24.data <- abs(DE_expression[,8] < 1.5) & abs(DE_expression[,10] < 1.5) & abs(DE_expression[,12] < 1.5) & DE_expression[,9] <0.05 & DE_expression[,11] <0.05 & DE_expression[,13] <0.05
table(WT18_24.vs.HET18_24.data)
M2_WT18_24.vs.HET18_24 <- M2[WT18_24.vs.HET18_24.data,]
rownames(M2_WT18_24.vs.HET18_24) <- BioM[WT18_24.vs.HET18_24.data,3]

# Data die hierboven gegenereerd is, wordt gezet in een heatmap.
# Deze heatmap bevat de DE data van verschillende comparisons.
colors <- colorpanel(50, "darkmagenta", "Cyan")
heatmap.2(M2_WT6.8_18_24.vs.HET6.8_18_24, ColSideColors = col_cell_age, density.info="none", col=colors, trace="none", scale = c("row"), cexRow=1)
heatmap.2(M2_WT18_24.vs.HET18_24, ColSideColors = col_cell_age, density.info="none", col=colors, trace="none", scale = c("row"), cexRow=1)

# Hieronder wordt alle data georderd en de eerste 50 genen met een hoge
# FDR expressie worden gebruikt om de expressie tussen de verschillende
# genen te laten zien.
M2.WT.HET <- M2[order(Toptable1[[1]]$FDR),][1:50,]
rownames(M2.WT.HET) <- BioM[order(Toptable1[[1]]$FDR),3][1:50]
M2.WT.HET.min6_8w <- M2[order(Toptable2[[1]]$FDR),][1:50,]
rownames(M2.WT.HET.min6_8w) <- BioM[order(Toptable1[[1]]$FDR),3][1:50]
M2.WT12.VS.HET18 <- M2[order(Toptable3[[1]]$FDR),][1:50,]
rownames(M2.WT12.VS.HET18) <- BioM[order(Toptable1[[1]]$FDR),3][1:50]

heatmap.2(M2.WT.HET, ColSideColors = col_cell_age, density.info="none", col=colors, trace="none", scale = c("row"), cexRow=0.6)
heatmap.2(M2.WT.HET.min6_8w, ColSideColors = col_cell_age, density.info="none", col=colors, trace="none", scale = c("row"), cexRow=0.6)
heatmap.2(M2.WT12.VS.HET18, ColSideColors = col_cell_age, density.info="none", col=colors, trace="none", scale = c("row"), cexRow=0.6)

####################################################################
#                          Venn Diagram                            #
####################################################################

vennDiagram(DE_expression[,c(15,17,19,3)] < 0.01 , cex=0.5, lwd=3)
earlyVisualization <- which(DE_expression[,15] < 0.01 & DE_expression[,17] < 0.01 & DE_expression[,19] < 0.01 & DE_expression[,3] < 0.01)
earlyVisualization <- DE_expression[earlyVisualization,]
sum(nrow(earlyVisualization))
earlyVisualization[,21]

vennDiagram(DE_expression[,c(5,7,9,3)] < 0.05, cex=0.5, lwd=3)
lateVisualization <- which(DE_expression[,5] < 0.05 & DE_expression[,7] < 0.05 & DE_expression[,9] < 0.05 & DE_expression[,3] < 0.05)
lateVisualization <- DE_expression[lateVisualization,]
sum(nrow(lateVisualization))
lateVisualization[,21]

vennDiagram(DE_expression[,c(11,13,3)] < 0.05, cex=0.5, circle.col = c(1:3), lwd=3)
WT.HETVisualization <- which(DE_expression[,11] < 0.05 & DE_expression[13] < 0.05 & DE_expression[,3] < 0.05)
WT.HETVisualization <- DE_expression[WT.HETVisualization,]
sum(nrow(WT.HETVisualization))
WT.HETVisualization[,21]
