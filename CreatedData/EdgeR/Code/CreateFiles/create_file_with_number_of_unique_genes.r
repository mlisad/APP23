####################################################################
# Author    : M. Dubbelaar
# Date      : 15-sept-2015
# File Name : create_file_with_number_of_unique_genes.r
# Purpose   : Creates a file with the number of unique genes of each DE.
####################################################################
# The unique gene names found within the toptable results after de FDR filtering are saved with the vector all.genes.

all.genes.linear.genotype <- read.table("../Made_Documents/All_ages/main_effect_genes_genotype.txt", sep = ",")
all.genes.linear.age <- read.table("../Made_Documents/All_ages/main_effect_genes_age.txt", sep=",")
all.genes.linear.interaction <- read.table("../Made_Documents/All_ages/linear_effect_genotype_age.txt", sep = ",")

olderMice.linear.genotype <- read.table("../Made_Documents/12-18-24M_old_mice/main_effect_genes_genotype_older_mice.txt", sep = ",")
olderMice.linear.age <- read.table("../Made_Documents/12-18-24M_old_mice/main_effect_genes_age_older_mice.txt", sep=",")
olderMice.linear.interaction <- read.table("../Made_Documents/12-18-24M_old_mice/linear_effect_genotype_age_older_mice.txt", sep=",")

youngerMice.linear.genotype <- read.table("../Made_Documents/6_8M-12M_old_mice/main_effect_genes_genotype_younger_mice.txt", sep=",")
youngerMice.linear.age <- read.table("../Made_Documents/6_8M-12M_old_mice/main_effect_genes_age_younger_mice.txt", sep=",")
youngerMice.linear.interaction <- read.table("../Made_Documents/6_8M-12M_old_mice/linear_effect_genotype_age_younger_mice.txt", sep=",")

all.genes <- NULL
all.genes <- unique(c(rownames(toptable1.results), rownames(toptable2.results), rownames(toptable3.results),
                      rownames(toptable4.results), rownames(toptable5.results), rownames(toptable6.results),
                      rownames(toptable7.results), rownames(toptable8.results), rownames(toptable9.results),
                      rownames(toptable10.results), rownames(toptable11.results), rownames(toptable12.results),
                      rownames(toptable13.results), rownames(toptable14.results), rownames(toptable15.results),
                      rownames(toptable16.results)))

# The number of the unique genes per comparison and of all genes are saved within a vector
unique.Genes2 <- rbind(length(rownames(toptable1.results)), length(rownames(toptable2.results)),
                       length(rownames(toptable3.results)), length(rownames(toptable4.results)),
                       length(rownames(toptable5.results)), length(rownames(toptable6.results)),
                       length(rownames(toptable7.results)), length(rownames(toptable8.results)),
                       length(rownames(toptable9.results)), length(rownames(toptable10.results)),
                       length(rownames(toptable11.results)), length(rownames(toptable12.results)),
                       length(rownames(toptable13.results)), length(rownames(toptable14.results)),
                       length(rownames(toptable15.results)), length(rownames(toptable16.results)),
                       length(all.genes), 
                       "", length(all.genes.linear.genotype[,1]), length(all.genes.linear.age[,1]), length(all.genes.linear.interaction[,1]), 
                       length(olderMice.linear.genotype[,1]), length(olderMice.linear.age[,1]), length(olderMice.linear.interaction[,1]), 
                       length(youngerMice.linear.genotype[,1]), length(youngerMice.linear.age[,1]), length(youngerMice.linear.interaction[,1]))

# The rownames will indicate with DE it is.
rownames(unique.Genes2) <-c("6.8W_WT-6.8W_HET", "12M_WT-12M_HET", "18M_WT-18M_HET", "24M_WT-24M_HET",
                            "12M_HET-6.8W_HET", "18M_HET-12M_HET", "24M_HET-18M_HET", "12M_WT-6.8W_WT",
                            "18M_WT-12M_WT", "24M_WT-18M_WT", "24M_HET-6.8W_HET", "24M_WT-6.8W_WT",
                            "(24M_WT-24M_HET)-(6.8W_WT-6.8W_HET)", "(12M_WT-12M_HET)-(6.8W_WT-6.8W_HET)",
                            "(18M_WT-18M_HET)-(12M_WT-12M_HET)", "(24M_WT-24M_HET)-(18M_WT-18M_HET)",
                            "All found genes",
                            "", "Main-effect genotype (ALL)", "Main-effect age (ALL)", "Interaction (ALL)",
                            "Main-effect genotype (olderMice)", "Main-effect age (olderMice)", "Interaction (olderMice)",
                            "Main-effect genotype (youngerMice)", "Main-effect age (youngerMice)", "Interaction (youngerMice")

# The data is written to a csv file, this can be used for further research.
write.table(unique.Genes2, "../Made_Documents/Unique_Genes_without_logFC.csv", col.names=F, sep=":")