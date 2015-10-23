####################################################################
# Author    : M. Dubbelaar
# Date      : 14-sept-2015
# File Name : Creation_main_MDS_plot.r 
# Purpose   : Creation different MDS plots with a color attached to the 
#             Phenotype.
####################################################################
####################################################################
#                      Plot for the age effect                     #
####################################################################

# A repetition will be made for the amount of dimentions of the target. 
col_cell_age <- rep("black", dim(targets)[1])
# The different conditions get a special color from dark to light.
col_cell_age
# Makes a red color for the WT mice.
# The colors go from light red (young mice) to dark red (old mice).
col_cell_age[targets$Conditie ==  "WT.06_08"] = col=rainbow(265)[20]
col_cell_age[targets$Conditie ==  "WT.12"] = col=rainbow(265)[15]
col_cell_age[targets$Conditie ==  "WT.18"] = col=rainbow(265)[10]
col_cell_age[targets$Conditie ==  "WT.24"] = col=rainbow(265)[1]
# Makes a blue color for the APP (defines as HET) mice.
# The colors go from light blue (young mice) to dark blue (old mice).
col_cell_age[targets$Conditie ==  "HET.06_08"] = col=rainbow(265)[145]
col_cell_age[targets$Conditie ==  "HET.12"] = col=rainbow(265)[150]
col_cell_age[targets$Conditie ==  "HET.18"] = col=rainbow(265)[155]
col_cell_age[targets$Conditie ==  "HET.24"] = col=rainbow(265)[160]
# A pdf will be made containing the plot, after this the creation will be stopped.
pdf("../Plots/MDS_plot_age_effect.pdf") 
?plotMDS.DGEList
plotMDS.DGEList(dge, col= col_cell_age, main="MDS plot (age effect)")
dev.off() 