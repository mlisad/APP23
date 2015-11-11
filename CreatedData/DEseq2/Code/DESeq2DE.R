####################################################################
# Author    : M. Dubbelaar
# Date      : 1-okt-2015
# File Name : DifferentialExpressionDESeq2.R
# Purpose   : To determine the unqiue genes within several DE.
####################################################################
# This condition needs to be available within the sampleTable, 
# to use this data in the design.
sampleCondition <- factor(targets$Conditie)
sampleTable <- data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
# DESeqDataSetFromHTSeqCount is used to store input values, intermediate .
# calculations and results of the DE analysis.
# The data will be filtered and normalized with the own made normalizeData function.
# The data will be used by the function DESeq after the normalisation.
# This function performs defauls analysis for the DE.
ddsHETSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~0+condition)
colData(ddsHETSeq)$condition <- factor(colData(ddsHETSeq)$condition, levels=unique(targets$Conditie))
dds <- normalizeData(ddsHETSeq)
dds <- DESeq(dds, betaPrior = F, full = ~0+condition)

# The different DE are shown below
# All of the conditions with an alpha of 0.05 will be saved.
WT2M.vs.HET2M <- results(dds, contrast = c("condition", "WT.02", "HET.02"), alpha = 0.05)
WT6M.vs.HET6M <- results(dds, contrast = c("condition", "WT.06", "HET.06"), alpha = 0.05)
WT18M.vs.HET18M <- results(dds, contrast = c("condition", "WT.18", "HET.18"), alpha = 0.05)
WT24M.vs.HET24M <- results(dds, contrast = c("condition", "WT.24", "HET.24"), alpha = 0.05)
HET6M.vs.HET2M <- results(dds, contrast = c("condition", "HET.06", "HET.02"), alpha = 0.05)
HET18M.vs.HET6M <- results(dds, contrast = c("condition", "HET.18", "HET.06"), alpha = 0.05)
HET24M.vs.HET18M <- results(dds, contrast = c("condition", "HET.24", "HET.18"), alpha = 0.05)
WT6M.vs.WT2M <- results(dds, contrast = c("condition", "WT.06","WT.02"), alpha = 0.05)
WT18M.vs.WT6M <- results(dds, contrast = c("condition", "WT.18", "WT.06"), alpha = 0.05)
WT24M.vs.WT18M <- results(dds, contrast = c("condition", "WT.24", "WT.18"), alpha = 0.05)
WT24M.vs.WT2M <- results(dds, contrast = c("condition","WT.24" , "WT.02"), alpha = 0.05)
HET24M.vs.HET2M <- results(dds, contrast = c("condition", "HET.24", "HET.02"), alpha = 0.05)

# The data will be used to get the unique genes with a adjusted p value of 0.05, 
# these genes will be saved.
toptable1.results <- getUniqueGenes(WT2M.vs.HET2M)
toptable2.results <- getUniqueGenes(WT6M.vs.HET6M)
toptable3.results <- getUniqueGenes(WT18M.vs.HET18M)
toptable4.results <- getUniqueGenes(WT24M.vs.HET24M)
toptable5.results <- getUniqueGenes(HET6M.vs.HET2M)
toptable6.results <- getUniqueGenes(HET18M.vs.HET6M)
toptable7.results <- getUniqueGenes(HET24M.vs.HET18M)
toptable8.results <- getUniqueGenes(WT6M.vs.WT2M)
toptable9.results <- getUniqueGenes(WT18M.vs.WT6M)
toptable10.results <- getUniqueGenes(WT24M.vs.WT18M)
toptable11.results <- getUniqueGenes(WT24M.vs.WT2M)
toptable12.results <- getUniqueGenes(HET24M.vs.HET2M)

