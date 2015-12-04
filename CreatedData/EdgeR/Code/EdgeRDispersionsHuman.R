####################################################################
# Author    : M. Dubbelaar
# Date      : 14-sept-2015
# File Name : edgeRDespersions.R
# Purpose   : Calculates the dispersions of each comparison
####################################################################
# The design is used with the calculation of the estimates. 
design <- model.matrix(~0+factor(targets$Genotype), data = dge$samples)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# A fit is made to read the counts for each transcript or gene.
fit <- glmFit(dge, design)

# Compares LOAD against WT.
lrt01 <- glmLRT(fit, contrast = c(-1,1))
####################################################################
#                           Filtering Data                         #
####################################################################
toptable1 <- createToptableResults(lrt01, "/home/mdubbelaar/Desktop/Results/Human/Made_Documents/DE/")
####################################################################
#                      Differential Expression                     #
####################################################################
DE.Expression <- cbind(rownames(toptable1[[1]][[1]]), toptable1[[1]][[1]]$logFC, toptable1[[1]][[1]]$FDR, BioM[,3:4])
geneColsHuman <- c("Genes",  "logFC: LOAD vs CTRL", "FDR: LOAD vs CTRL", "Gene Symbol", "Gene Description")
write.table(DE.Expression , "/home/mdubbelaar/Desktop/Results/Human/Made_Documents/DE_Files/DifferentialExpressionGenotypes", row.names = F, col.names = geneColsHuman, sep = "\t")
