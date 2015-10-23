####################################################################
# Author    : M. Dubbelaar
# Date      : 2-okt-2015
# File Name : DESeQ2_uniqueGenesFile.r
# Purpose   : Creates a txt file with the amount of unique genes for each DE.
####################################################################
all.genes.linear.genotype <- read.table("Made_Documents/All_ages/main_genotype_results.txt", sep = ",")
all.genes.linear.age <- read.table("Made_Documents/All_ages/main_age_results.txt", sep=",")
all.genes.linear.interaction <- read.table("Made_Documents/All_ages/interaction_results.txt", sep = ",")

olderMice.linear.genotype <- read.table("Made_Documents/12-18-24M_old_mice/main_genotype_results.txt", sep = ",")
olderMice.linear.age <- read.table("Made_Documents/12-18-24M_old_mice/main_age_results.txt", sep=",")
olderMice.linear.interaction <- read.table("Made_Documents/12-18-24M_old_mice/interaction_results.txt", sep=",")

youngerMice.linear.genotype <- read.table("Made_Documents/6_8M-12M_old_mice/main_genotype_results.txt", sep=",")
youngerMice.linear.age <- read.table("Made_Documents/6_8M-12M_old_mice/main_age_results.txt", sep=",")
youngerMice.linear.interaction <- read.table("Made_Documents/6_8M-12M_old_mice/interaction_results.txt", sep=",")

all.genes <- NULL
all.genes <- unique(c(rownames(WT06_08.vs.HET06_08), rownames(WT12M.vs.HET12M), rownames(WT18M.vs.HET18M),
                      rownames(WT24M.vs.HET24M), rownames(HET12M.vs.HET06_08M), rownames(HET18M.vs.HET12M),
                      rownames(HET24M.vs.HET18M), rownames(WT12M.vs.WT06_08M), rownames(WT18M.vs.WT12M),
                      rownames(WT24M.vs.WT18M), rownames(WT24M.vs.WT06_08W), rownames(HET24M.vs.HET06_08W)))

# The number of the unique genes per comparison and of all genes are saved within a vector
unique.Genes2 <- rbind(length(rownames(WT06_08.vs.HET06_08)), length(rownames(WT12M.vs.HET12M)),
                       length(rownames(WT18M.vs.HET18M)), length(rownames(WT24M.vs.HET24M)),
                       length(rownames(HET12M.vs.HET06_08M)), length(rownames(HET18M.vs.HET12M)),
                       length(rownames(HET24M.vs.HET18M)), length(rownames(WT12M.vs.WT06_08M)),
                       length(rownames(WT18M.vs.WT12M)), length(rownames(WT24M.vs.WT18M)),
                       length(rownames(WT24M.vs.WT06_08W)), length(rownames(HET24M.vs.HET06_08W)),
                       length(all.genes), 
                       "", length(all.genes.linear.genotype[,1]), length(all.genes.linear.age[,1]), length(all.genes.linear.interaction[,1]), 
                       length(olderMice.linear.genotype[,1]), length(olderMice.linear.age[,1]), length(olderMice.linear.interaction[,1]), 
                       length(youngerMice.linear.genotype[,1]), length(youngerMice.linear.age[,1]), length(youngerMice.linear.interaction[,1]))

# The rownames will indicate with DE it is.
rownames(unique.Genes2) <-c("6.8W_WT-6.8W_HET", "12M_WT-12M_HET", "18M_WT-18M_HET", "24M_WT-24M_HET",
                            "12M_HET-6.8W_HET", "18M_HET-12M_HET", "24M_HET-18M_HET", "12M_WT-6.8W_WT",
                            "18M_WT-12M_WT", "24M_WT-18M_WT", "24M_WT-6.8W_WT", "24M_HET-6.8W_HET", 
                            "All found genes",
                            "", "Main-effect genotype (ALL)", "Main-effect age (ALL)", "Interaction (ALL)",
                            "Main-effect genotype (olderMice)", "Main-effect age (olderMice)", "Interaction (olderMice)",
                            "Main-effect genotype (youngerMice)", "Main-effect age (youngerMice)", "Interaction (youngerMice)")

# The data is written to a csv file, this can be used for further research.
write.table(unique.Genes2, "Made_Documents/Unique_Genes_with_0.05_AdjpVal.csv", col.names=F, sep=":")
