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
# Compares 6 month WT against 6 month HET
lrt02 <- glmLRT(fit, contrast = c(0,1,0,0,0,-1,0,0))
# Compares 18 month WT against 18 month HET
lrt03 <- glmLRT(fit, contrast = c(0,0,1,0,0,0,-1,0))
# Compares 24 month WT against 24 month HET
lrt04 <- glmLRT(fit, contrast = c(0,0,0,1,0,0,0,-1))

# Compares 6 month HET against 6-8 week HET
lrt05 <- glmLRT(fit, contrast = c(1,-1,0,0,0,0,0,0))
# Compares 18 month HET against 6 month HET
lrt06 <- glmLRT(fit, contrast = c(0,1,-1,0,0,0,0,0))
# Compares 24 month HET against 18 month HET
lrt07 <- glmLRT(fit, contrast = c(0,0,1,-1,0,0,0,0))

# Compares 6 month  WT against 6-8 week WT
lrt08 <- glmLRT(fit, contrast = c(0,0,0,0,1,-1,0,0))
# Compares 18 month WT against 6 month WT
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
#                            Heatmaps                              #
####################################################################
# Checks the two toptable results with the most unique genes.
# These heatmaps are saved within the heatmap.pdf.
#pdf("/home/mdubbelaar/Desktop/Results/APP23/Plots/heatmaps_containing_most_genes.pdf") 
pdf("/Volumes/Elements_Marissa/School/Stage/APP23/Results/APP23_results/EdgeR/Plots/heatmaps_containing_most_genes.pdf") 
heatmap.2(M2[match(rownames(toptable11.results), rownames(M2)),c(4:6, 10:12, 16:18, 22:24)], ColSideColors = col_cell_age[c(4:6, 10:12, 16:18, 22:24)], cexRow = 0.01, trace = "none", scale = "row", main="24M HET - 2M HET")
heatmap.2(M2[match(rownames(toptable12.results), rownames(M2)),c(1:3, 7:9, 13:15, 19:21)], ColSideColors = col_cell_age[c(1:3, 7:9, 13:15, 19:21)], cexRow = 0.01, trace = "none", scale = "row", main="24M WT - 2M WT")
dev.off()

####################################################################
#                      Differential Expression                     #
####################################################################
DE.ExpressionWT <- cbind(rownames(toptable1[[1]][[1]]), toptable8[[1]][[1]]$logFC, toptable8[[1]][[1]]$FDR, toptable9[[1]][[1]]$logFC, toptable9[[1]][[1]]$FDR,
                         toptable10[[1]][[1]]$logFC, toptable10[[1]][[1]]$FDR, toptable12[[1]][[1]]$logFC, toptable12[[1]][[1]]$FDR, BioM[,3:4])

DE.ExpressionHET <- cbind(rownames(toptable1[[1]][[1]]),  toptable5[[1]][[1]]$logFC, toptable5[[1]][[1]]$FDR, toptable6[[1]][[1]]$logFC, toptable6[[1]][[1]]$FDR,
                         toptable7[[1]][[1]]$logFC, toptable7[[1]][[1]]$FDR, toptable11[[1]][[1]]$logFC, toptable11[[1]][[1]]$FDR, BioM[,3:4])
  
DE.ExpressionCombi <- cbind(rownames(toptable1[[1]][[1]]), toptable1[[1]][[1]]$logFC, toptable1[[1]][[1]]$FDR, toptable2[[1]][[1]]$logFC, toptable2[[1]][[1]]$FDR,
                            toptable3[[1]][[1]]$logFC, toptable3[[1]][[1]]$FDR, toptable4[[1]][[1]]$logFC, toptable4[[1]][[1]]$FDR, BioM[,3:4])
####################################################################
#                         Creating DE Files                        #
####################################################################
geneColsWT <- c("Genes",  "logFC: 6M WT vs 2M WT", "FDR: 6M WT vs 2M WT", "logFC: 18M WT vs 6M WT", "FDR: 18M WT vs 6M WT", 
                "logFC: 24M WT vs 18M WT", "FDR: 24M WT vs 18M WT", "logFC: 24M WT vs 2M WT", "FDR: 24M WT vs 2M WT", "Gene Symbol", "Gene Description")
geneColsHET <- c("Genes", "logFC: 6M HET vs 2M HET", "FDR: 6M HET vs 2M HET", "logFC: 18M HET vs 6M HET", "FDR: 18M HET vs 6M HET", 
                 "logFC: 24M HET vs 18M HET", "FDR: 24M HET vs 18M HET", "logFC: 24M HET vs 2M HET", "FDR: 24M HET vs 2M HET", "Gene Symbol", "Gene Description")
geneColsCombi <- c("Genes", "logFC: 2M WT vs 2M HET","FDR: 2M WT vs 2M HET", "logFC: 6M WT vs 6M HET",  "FDR: 6M WT vs 6M HET", 
                   "logFC: 18M WT vs 18M HET", "FDR: 18M WT vs 18M HET", "logFC: 24M WT vs 24M HET", "FDR: 24M WT vs 24M HET", "Gene Symbol", "Gene Description")
#write.table(DE.ExpressionWT , "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesWT.txt", row.names = F, col.names = geneColsWT, sep = "\t")
#write.table(DE.ExpressionHET , "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesHET.txt", row.names = F, col.names = geneColsHET, sep = "\t")
#write.table(DE.ExpressionCombi , "/home/mdubbelaar/Desktop/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesHET-WT.txt", row.names = F, col.names = geneColsCombi, sep = "\t")
write.table(DE.ExpressionWT , "/Volumes/Elements_Marissa/School/Stage/APP23/Results/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesWT.txt", row.names = F, col.names = geneColsWT, sep = "\t")
write.table(DE.ExpressionHET , "/Volumes/Elements_Marissa/School/Stage/APP23/Results/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesHET.txt", row.names = F, col.names = geneColsHET, sep = "\t")
write.table(DE.ExpressionCombi , "/Volumes/Elements_Marissa/School/Stage/APP23/Results/APP23_results/EdgeR/Made_Documents/DE_Files/DifferentialGenesHET-WT.txt", row.names = F, col.names = geneColsCombi, sep = "\t")

