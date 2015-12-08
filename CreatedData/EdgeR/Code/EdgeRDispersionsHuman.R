####################################################################
# Author    : M. Dubbelaar
# Date      : 05-dec-2015
# File Name : edgeRDespersionsHuman.R
# Purpose   : Calculates the differential expression genes among
#             the LOAD and the CTRL patients.
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
#                         Creating DE Files                        #
####################################################################
DE.Expression <- cbind(rownames(toptable1[[1]][[1]]), toptable1[[1]][[1]]$logFC, toptable1[[1]][[1]]$FDR, BioM[,3:4])
geneColsHuman <- c("Genes",  "logFC: LOAD vs CTRL", "FDR: LOAD vs CTRL", "Gene Symbol", "Gene Description")
write.table(DE.Expression , "/home/mdubbelaar/Desktop/Results/Human/Made_Documents/DE_Files/DifferentialExpressionGenotypes", row.names = F, col.names = geneColsHuman, sep = "\t")
