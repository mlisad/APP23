####################################################################
# Author    : M. Dubbelaar
# Date      : 04-dec-2015
# File Name : edgeRDespersionsCKp25.R
# Purpose   : Calculates the dispersions of each comparison
####################################################################
# The design is used with the calculation of the estimates. 
design <- model.matrix(~0+factor(targets$Conditie), data = dge$samples)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)



# A fit is made to read the counts for each transcript or gene.
fit <- glmFit(dge, design)

# Compares 2 week WT against 6 week WT
lrt01 <- glmLRT(fit, contrast = c(0,0,1,-1))
# Compares 2 week WT against 2 week CKp25
lrt02 <- glmLRT(fit, contrast = c(-1,0,1,0))
# Compares 2 week WT against 6 week CKp25
lrt03 <- glmLRT(fit, contrast = c(0,-1,1,0))
# Compares 6 week WT against 2 week CKp25
lrt04 <- glmLRT(fit, contrast = c(-1,0,0,1))
# Compares 6 week WT against 6 week CKp25
lrt05 <- glmLRT(fit, contrast = c(0,-1,0,1))
# Compares 2 week CKp25 against 6 week CKp25
lrt06 <- glmLRT(fit, contrast = c(1,-1,0,0))
####################################################################
#                           Filtering Data                         #
####################################################################
toptable1 <- createToptableResults(lrt01, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
toptable2 <- createToptableResults(lrt02, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
toptable3 <- createToptableResults(lrt03, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
toptable4 <- createToptableResults(lrt04, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
toptable5 <- createToptableResults(lrt05, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
toptable6 <- createToptableResults(lrt06, "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE/")
####################################################################
#                      Differential Expression                     #
####################################################################

DE.Expression <- cbind(rownames(toptable1[[1]][[1]]), toptable1[[1]][[1]]$logFC, toptable1[[1]][[1]]$FDR, toptable2[[1]][[1]]$logFC, toptable2[[1]][[1]]$FDR,
                         toptable3[[1]][[1]]$logFC, toptable3[[1]][[1]]$FDR, toptable4[[1]][[1]]$logFC, toptable4[[1]][[1]]$FDR,
                         toptable5[[1]][[1]]$logFC, toptable5[[1]][[1]]$FDR, toptable6[[1]][[1]]$logFC, toptable6[[1]][[1]]$FDR, BioM[,3:4])
geneCols <- c("Genes",  "logFC: 2W WT vs 6W WT", "FDR: 2W WT vs 6W WT", "logFC: 2W WT vs 2W CKp25", "FDR: 2W WT vs 2W CKp25", 
                "logFC: 2W WT vs 6W CKp25", "FDR: 2W WT vs 6W CKp25", "logFC: 6W WT vs 2W CKp25", "FDR: 6W WT vs 2W CKp25", 
                "logFC: 6W WT vs 6W CKp25", "FDR: 6W WT vs 6W CKp25", "logFC: 2W CKp25 vs 6W CKp25", "FDR: 2W CKp25 vs 6W CKp25","Gene Symbol", "Gene Description")
write.table(DE.Expression , "/home/mdubbelaar/Desktop/Results/CKp25/Made_Documents/DE_Files/DifferentialGenesCKp25.txt", row.names = F, col.names = geneCols, sep = "\t")

