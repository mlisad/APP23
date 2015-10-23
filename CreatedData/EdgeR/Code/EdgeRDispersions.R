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

# With the use of decideTestsDGE a detection is made to distinguish the up, down al all regulated genes.
table(decideTestsDGE(lrt01, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt02, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt03, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt04, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt05, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt06, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt07, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt08, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt09, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt10, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt11, p=0.05, adjust="BH"))
table(decideTestsDGE(lrt12, p=0.05, adjust="BH"))

# Data is stored within the toptables for further anaylsis.
toptable1 <- topTags(lrt01, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable2 <- topTags(lrt02, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable3 <- topTags(lrt03, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable4 <- topTags(lrt04, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable5 <- topTags(lrt05, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable6 <- topTags(lrt06, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable7 <- topTags(lrt07, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable8 <- topTags(lrt08, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable9 <- topTags(lrt09, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable10 <- topTags(lrt10, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable11 <- topTags(lrt11, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")
toptable12 <- topTags(lrt12, n=dim(dge[[1]])[1], adjust.method="BH", sort.by="none")

# The results from the toptables are checked with the use of the own-made filtergenes function
# This function returns the data which contain a FDR value below the 0.05
toptable1.results <- filterGenes(toptable1)
toptable2.results <- filterGenes(toptable2)
toptable3.results <- filterGenes(toptable3)
toptable4.results <- filterGenes(toptable4)
toptable5.results <- filterGenes(toptable5)
toptable6.results <- filterGenes(toptable6)
toptable7.results <- filterGenes(toptable7)
toptable8.results <- filterGenes(toptable8)
toptable9.results <- filterGenes(toptable9)
toptable10.results <- filterGenes(toptable10)
toptable11.results <- filterGenes(toptable11)
toptable12.results <- filterGenes(toptable12)