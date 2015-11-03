####################################################################
# Author    : M. Dubbelaar
# Date      : 29-sept-2015
# File Name : saveDEData.R
# Purpose   : Creates a txt file with the amount of unique genes 
#             for the DE.
####################################################################
# The data of the effects found with the the linear time design
# need to be read from the saved files. These files are the younger
# mice, the older mice and mice of all ages.


allGenesLinearInteraction <- NULL
youngerMiceLinearInteraction <- NULL

allGenesLinearGenotype <- read.table("Made_Documents/All_ages/main_genotype_result.txt", sep = ",")
allGenesLinearAge <- read.table("Made_Documents/All_ages/main_age_result.txt", sep=",")
if (length(readLines("Made_Documents/All_ages/interaction_result.txt"))-1 == 0) {
  allGenesLinearInteraction <- read.table("Made_Documents/All_ages/interaction_result.txt", sep = ",")
} 
#allGenesLinearInteraction <- read.table("Made_Documents/All_ages/interaction_result.txt", sep = ",")
olderMiceLinearGenotype <- read.table("Made_Documents/12-18-24M_old_mice/main_genotype_result.txt", sep = ",")
olderMiceLinearAge <- read.table("Made_Documents/12-18-24M_old_mice/main_age_result.txt", sep=",")
olderMiceLinearInteraction <- read.table("Made_Documents/12-18-24M_old_mice/interaction_result.txt", sep=",")
youngerMiceLinearGenotype <- read.table("Made_Documents/6_8M-12M_old_mice/main_genotype_result.txt", sep=",")
youngerMiceLinearAge <- read.table("Made_Documents/6_8M-12M_old_mice/main_age_result.txt", sep=",")
if (length(readLines("Made_Documents/6_8M-12M_old_mice/interaction_result.txt"))-1 ==0) {
  youngerMiceLinearInteraction <- read.table("Made_Documents/6_8M-12M_old_mice/interaction_result.txt", sep=",")
} 

# The unique genenames found within the normal DE are saved with 
# the vector allGenes.
allGenes <- NULL
allGenes <- unique(c(rownames(toptable1.results), rownames(toptable2.results), rownames(toptable3.results),
                      rownames(toptable4.results), rownames(toptable5.results), rownames(toptable6.results),
                      rownames(toptable7.results), rownames(toptable8.results), rownames(toptable9.results),
                      rownames(toptable10.results), rownames(toptable11.results), rownames(toptable12.results)))
####################################################################
# One of the packages lack genes within the youngerMiceLinearInteraction,
# this is the reason of the if - else usage. The length of the known 
# vectors are used to create a matrix with the number of unique genes
# for each table.
uniqueGenes <- NULL
if (!is.null(youngerMiceLinearInteraction) & !is.null(allGenesLinearInteraction)) {
  unique.Genes <- rbind(length(rownames(toptable1.results)), length(rownames(toptable2.results)),
                         length(rownames(toptable3.results)), length(rownames(toptable4.results)),
                         length(rownames(toptable5.results)), length(rownames(toptable6.results)),
                         length(rownames(toptable7.results)), length(rownames(toptable8.results)),
                         length(rownames(toptable9.results)), length(rownames(toptable10.results)),
                         length(rownames(toptable11.results)), length(rownames(toptable12.results)), length(allGenes), 
                         "", length(allGenesLinearGenotype[,1]), length(allGenesLinearAge[,1]), length(allGenesLinearInteraction[,1]), 
                         length(olderMiceLinearGenotype[,1]), length(olderMiceLinearAge[,1]), length(olderMiceLinearInteraction[,1]), 
                         length(youngerMiceLinearGenotype[,1]), length(youngerMiceLinearAge[,1]), length(youngerMiceLinearInteraction[,1]))

# The rownames will indicate with DE it is.
rownames(unique.Genes) <-c("6.8W_WT-6.8W_HET", "12M_WT-12M_HET", "18M_WT-18M_HET", "24M_WT-24M_HET",
                            "12M_HET-6.8W_HET", "18M_HET-12M_HET", "24M_HET-18M_HET", "12M_WT-6.8W_WT",
                            "18M_WT-12M_WT", "24M_WT-18M_WT", "24M_HET-6.8W_HET", "24M_WT-6.8W_WT", "All found genes",
                            "", "Main-effect genotype (ALL)", "Main-effect age (ALL)", "Interaction (ALL)",
                            "Main-effect genotype (olderMice)", "Main-effect age (olderMice)", "Interaction (olderMice)",
                            "Main-effect genotype (youngerMice)", "Main-effect age (youngerMice)", "Interaction (youngerMice)")
} else {
  unique.Genes <- rbind(length(rownames(toptable1.results)), length(rownames(toptable2.results)),
                        length(rownames(toptable3.results)), length(rownames(toptable4.results)),
                        length(rownames(toptable5.results)), length(rownames(toptable6.results)),
                        length(rownames(toptable7.results)), length(rownames(toptable8.results)),
                        length(rownames(toptable9.results)), length(rownames(toptable10.results)),
                        length(rownames(toptable11.results)), length(rownames(toptable12.results)), length(allGenes), 
                        "", length(allGenesLinearGenotype[,1]), length(allGenesLinearAge[,1]), 
                        length(olderMiceLinearGenotype[,1]), length(olderMiceLinearAge[,1]), length(olderMiceLinearInteraction[,1]), 
                        length(youngerMiceLinearGenotype[,1]), length(youngerMiceLinearAge[,1]))
  
  rownames(unique.Genes) <-c("6.8W_WT-6.8W_HET", "12M_WT-12M_HET", "18M_WT-18M_HET", "24M_WT-24M_HET",
                              "12M_HET-6.8W_HET", "18M_HET-12M_HET", "24M_HET-18M_HET", "12M_WT-6.8W_WT",
                              "18M_WT-12M_WT", "24M_WT-18M_WT", "24M_HET-6.8W_HET", "24M_WT-6.8W_WT", "All found genes",
                              "", "Main-effect genotype (ALL)", "Main-effect age (ALL)",
                              "Main-effect genotype (olderMice)", "Main-effect age (olderMice)", "Interaction (olderMice)",
                              "Main-effect genotype (youngerMice)", "Main-effect age (youngerMice)")
}
####################################################################
# The data is written to a csv file, this can be used for further research.
write.table(unique.Genes, "Made_Documents/Amount_Differential_Expression.csv", col.names=F, sep=":")