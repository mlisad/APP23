####################################################################
# Author    : M. Dubbelaar
# Date      : 14-sept-2015
# File Name : edgeRDespersions.R
# Purpose   : Calculates the dispersions of each comparison
####################################################################
# The design is used with the calculation of the estimates. 
design <- model.matrix(~0+factor(targets$Conditie), data = dge$samples)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

# A fit is made to read the counts for each transcript or gene.
fit <- glmFit(dge, design)

# Compares 6-8 week WT against 6-8 week HET
lrt01 <- glmLRT(fit, contrast = c(1,0,0,0,-1,0,0,0))
# Compares 12 month WT against 12 month HET
lrt02 <- glmLRT(fit, contrast = c(0,1,0,0,0,-1,0,0))
# Compares 18 month WT against 18 month HET
lrt03 <- glmLRT(fit, contrast = c(0,0,1,0,0,0,-1,0))
# Compares 24 month WT against 24 month HET
lrt04 <- glmLRT(fit, contrast = c(0,0,0,1,0,0,0,-1))

# Compares 12 month HET against 6-8 week HET
lrt05 <- glmLRT(fit, contrast = c(1,-1,0,0,0,0,0,0))
# Compares 18 month HET against 12 month HET
lrt06 <- glmLRT(fit, contrast = c(0,1,-1,0,0,0,0,0))
# Compares 24 month HET against 18 month HET
lrt07 <- glmLRT(fit, contrast = c(0,0,1,-1,0,0,0,0))

# Compares 12 month  WT against 6-8 week WT
lrt08 <- glmLRT(fit, contrast = c(0,0,0,0,1,-1,0,0))
# Compares 18 month WT against 12 month WT
lrt09 <- glmLRT(fit, contrast = c(0,0,0,0,0,1,-1,0))
# Compares 24 month WT against 18 month WT
lrt10 <- glmLRT(fit, contrast = c(0,0,0,0,0,0,1,-1))

# Compares 24 month HET against 6-8 week HET
lrt11 <- glmLRT(fit, contrast = c(1,0,0,-1,0,0,0,0))
# Compares 24 month WT against 6-8 week WT
lrt12 <- glmLRT(fit, contrast = c(0,0,0,0,1,0,0,-1))

####################################################################
#                           Filtering Data                         #
####################################################################

createToptableResults <- function(lrt, number){
  # With the use of decideTestsDGE a detection is made to distinguish the up, down al all regulated genes.
  print(table(decideTestsDGE(lrt, p=0.05, adjust="BH")))
  # Data is stored within the toptables for further anaylsis.
  toptable <- topTags(lrt, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  # The results from the toptables are checked with the use of the own-made filtergenes function
  # This function returns the data which contain a FDR value below the 0.05
  toptable.results <- filterGenes(toptable)
  #write.table(rownames(toptable.results), paste("Made_Documents/DE/DEtable", number, sep=""), sep="\n", row.names = F, col.names = F)
  return (toptable.results)
} 

toptable1.results <- createToptableResults(lrt01, "1")
toptable2.results <- createToptableResults(lrt02, "2")
toptable3.results <- createToptableResults(lrt03, "3")
toptable4.results <- createToptableResults(lrt04, "4")
toptable5.results <- createToptableResults(lrt05, "5")
toptable6.results <- createToptableResults(lrt06, "6")
toptable7.results <- createToptableResults(lrt07, "7")
toptable8.results <- createToptableResults(lrt08, "8")
toptable9.results <- createToptableResults(lrt09, "9")
toptable10.results <- createToptableResults(lrt10, "10")
toptable11.results <- createToptableResults(lrt11, "11")
toptable12.results <- createToptableResults(lrt12, "12")
####################################################################
#                      Differential Expression                     #
####################################################################
DE.ExpressionLogFC <- cbind(rownames(toptable1[[1]]), toptable1[[1]]$logFC, toptable2[[1]]$logFC, toptable3[[1]]$logFC,
                       toptable4[[1]]$logFC, toptable5[[1]]$logFC, toptable6[[1]]$logFC,
                       toptable7[[1]]$logFC, toptable8[[1]]$logFC, toptable9[[1]]$logFC,
                       toptable10[[1]]$logFC, toptable11[[1]]$logFC, toptable12[[1]]$logFC, BioM[,3:4])
DE.ExpressionFDR <- cbind(rownames(toptable1[[1]]), toptable1[[1]]$FDR, toptable2[[1]]$FDR, toptable3[[1]]$FDR,
                            toptable4[[1]]$FDR, toptable5[[1]]$FDR, toptable6[[1]]$FDR,
                            toptable7[[1]]$FDR, toptable8[[1]]$FDR, toptable9[[1]]$FDR,
                            toptable10[[1]]$FDR, toptable11[[1]]$FDR, toptable12[[1]]$FDR, BioM[,3:4])

####################################################################
#                         Creating DE Files                        #
####################################################################
geneCols <- c("Genes", "2M WT vs 2M HET", "12M WT vs 12M HET", "18M WT vs 18M HET", "24M WT vs 24M HET", 
              "12M HET vs 2M HET", "18M HET vs 12M HET", "24M HET vs 18M HET", "12M WT vs 2M WT", "18M WT vs 12M WT", "24M WT vs 18M WT",
              "24M HET vs 2M HET", "24M WT vs 2M WT","Gene Symbol", "Gene Description")
write.table(DE.ExpressionLogFC , "Made_Documents/DifferentialGenesLogFC.txt", row.names = F, col.names = geneCols, sep = "\t")
write.table(DE.ExpressionFDR , "Made_Documents/DifferentialGenesFDR.txt", row.names = F, col.names = geneCols, sep = "\t")

