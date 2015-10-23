####################################################################
# Author    : M. Dubbelaar
# Date      : 18-sept-2015
# File Name : meer_testen.r
# Purpose   : 
# Used Files: Visualization_APP23_data.R
####################################################################

# Nc zijn de genormaliseerde columns
nc <- cpm(dge, normalized.lib.sizes = T)
# rn zijn alleen de rownames van de genen in toptable.
rn1 <- rownames(Toptable1$table)
rn2 <- rownames(Toptable2$table)
rn3 <- rownames(Toptable3$table)
rn4 <- rownames(Toptable4$table)
rn5 <- rownames(Toptable5$table)
rn6 <- rownames(Toptable6$table)
rn7 <- rownames(Toptable7$table)
rn8 <- rownames(Toptable8$table)
rn9 <- rownames(Toptable9$table)
# Hieronder zijn de genen te vinden met de cpm waarden.
# De columns zijn gerangschikt op de conditie van de samples.
head(nc[rn1,order(targets$Conditie)], 5)
head(nc[rn2,order(targets$Conditie)], 5)
head(nc[rn3,order(targets$Conditie)], 5)
head(nc[rn4,order(targets$Conditie)], 5)
head(nc[rn5,order(targets$Conditie)], 5)
head(nc[rn6,order(targets$Conditie)], 5)
head(nc[rn7,order(targets$Conditie)], 5)
head(nc[rn8,order(targets$Conditie)], 5)
head(nc[rn9,order(targets$Conditie)], 5)

deg1 <- rn1[Toptable1$table$FDR < .05]
deg2 <- rn2[Toptable2$table$FDR < .05]
deg3 <- rn3[Toptable3$table$FDR < .05]
deg4 <- rn4[Toptable4$table$FDR < .05]
deg5 <- rn5[Toptable5$table$FDR < .05]
deg6 <- rn6[Toptable6$table$FDR < .05]
deg7 <- rn7[Toptable7$table$FDR < .05]
deg8 <- rn8[Toptable8$table$FDR < .05]
deg9 <- rn9[Toptable9$table$FDR < .05]

# Mbv de plotsmear kan er worden gekeken welke genen (met een p value van 5%) tot expressie komen.
plotSmear(dge, de.tags = deg1)
plotSmear(dge, de.tags = deg2)
plotSmear(dge, de.tags = deg3)
plotSmear(dge, de.tags = deg4)
plotSmear(dge, de.tags = deg5)
plotSmear(dge, de.tags = deg6)
plotSmear(dge, de.tags = deg7)
plotSmear(dge, de.tags = deg8)
plotSmear(dge, de.tags = deg9)



