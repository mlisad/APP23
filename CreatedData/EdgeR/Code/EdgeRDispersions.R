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
createToptableResults <- function(lrt){
  # With the use of decideTestsDGE a detection is made to distinguish the up, down al all regulated genes.
  print(table(decideTestsDGE(lrt, p=0.05, adjust="BH")))
  # Data is stored within the toptables for further anaylsis.
  toptable <- topTags(lrt, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
  # The results from the toptables are checked with the use of the own-made filtergenes function
  # This function returns the data which contain a FDR value below the 0.05
  toptable.results <- filterGenes(toptable)
  # The comparison of the samples will be split at the ")" sign for the first time,
  # the second time the 2th string is split at the " -".
  # this results two genotype and age sets.
  # These information will be added together to create the filename.
  firstSplit <- strsplit(lrt$comparison, split = ")")
  secondSplit <- strsplit(firstSplit[[1]][2], split = " -")
  write.table(rownames(toptable.results), paste("Made_Documents/DE/DE_", firstSplit[[1]][3], "_vs_", secondSplit[[1]][1], ".txt", sep = ""), eol=",\n", quote = F, row.names = F, col.names = F)
  output <- list(toptable, toptable.results)
  return (output)
} 

toptable1 <- createToptableResults(lrt01)
toptable1.results <- toptable1[[2]]
toptable2 <- createToptableResults(lrt02)
toptable2.results <- toptable2[[2]]
toptable3 <- createToptableResults(lrt03)
toptable3.results <- toptable3[[2]]
toptable4 <- createToptableResults(lrt04)
toptable4.results <- toptable4[[2]]
toptable5 <- createToptableResults(lrt05)
toptable5.results <- toptable5[[2]]
toptable6 <- createToptableResults(lrt06)
toptable6.results <- toptable6[[2]]
toptable7 <- createToptableResults(lrt07)
toptable7.results <- toptable7[[2]]
toptable8 <- createToptableResults(lrt08)
toptable8.results <- toptable8[[2]]
toptable9 <- createToptableResults(lrt09)
toptable9.results <- toptable9[[2]]
toptable10 <- createToptableResults(lrt10)
toptable10.results <- toptable10[[2]]
toptable11 <- createToptableResults(lrt11)
toptable11.results <- toptable11[[2]]
toptable12 <- createToptableResults(lrt12)
toptable12.results <- toptable12[[2]]
####################################################################
#                      Differential Expression                     #
####################################################################

DE.ExpressionLogFC <- cbind(rownames(toptable1.results[[1]][[1]]), toptable1.results[[1]][[1]]$logFC, toptable2.results[[1]][[1]]$logFC, toptable3.results[[1]][[1]]$logFC,
                            toptable4.results[[1]][[1]]$logFC, toptable5.results[[1]][[1]]$logFC, toptable6.results[[1]][[1]]$logFC,
                            toptable7.results[[1]][[1]]$logFC, toptable8.results[[1]][[1]]$logFC, toptable9.results[[1]][[1]]$logFC,
                            toptable10.results[[1]][[1]]$logFC, toptable11.results[[1]][[1]]$logFC, toptable12.results[[1]][[1]]$logFC, BioM[,3:4])
DE.ExpressionFDR <- cbind(rownames(toptable1.results[[1]][[1]]), toptable1.results[[1]][[1]]$FDR, toptable2.results[[1]][[1]]$FDR, toptable3.results[[1]][[1]]$FDR,
                          toptable4.results[[1]][[1]]$FDR, toptable5.results[[1]][[1]]$FDR, toptable6.results[[1]][[1]]$FDR,
                          toptable7.results[[1]][[1]]$FDR, toptable8.results[[1]][[1]]$FDR, toptable9.results[[1]][[1]]$FDR,
                          toptable10.results[[1]][[1]]$FDR, toptable11.results[[1]][[1]]$FDR, toptable12.results[[1]][[1]]$FDR, BioM[,3:4])
####################################################################
#                         Creating DE Files                        #
####################################################################
geneCols <- c("Genes", "2M WT vs 2M HET", "12M WT vs 12M HET", "18M WT vs 18M HET", "24M WT vs 24M HET", 
              "12M HET vs 2M HET", "18M HET vs 12M HET", "24M HET vs 18M HET", "12M WT vs 2M WT", "18M WT vs 12M WT", "24M WT vs 18M WT",
              "24M HET vs 2M HET", "24M WT vs 2M WT", "Gene Symbol", "Gene Description")
write.table(DE.ExpressionLogFC , "Made_Documents/DifferentialGenesLogFC.txt", row.names = F, col.names = geneCols, sep = "\t")
write.table(DE.ExpressionFDR , "Made_Documents/DifferentialGenesFDR.txt", row.names = F, col.names = geneCols, sep = "\t")

