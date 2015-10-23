####################################################################
# Author    : M. Dubbelaar
# Date      : 14-sept-2015
# File Name : Creation_MDS_plots.r 
# Purpose   : Creation different MDS plots with a color attached to the 
#             Phenotype.
####################################################################
####################################################################
#                        Plot for the plate                        #
####################################################################
col_cell_plate <- rep("black", dim(targets)[1])
col_cell_plate[targets$Plaat.nr ==  "A1"] = "navy"
col_cell_plate[targets$Plaat.nr ==  "A2"] = "blue3"
col_cell_plate[targets$Plaat.nr ==  "A3"] = "mediumblue"
col_cell_plate[targets$Plaat.nr ==  "A4"] = "blue1"
col_cell_plate[targets$Plaat.nr ==  "A5"] = "dodgerblue"
col_cell_plate[targets$Plaat.nr ==  "A6"] = "deepskyblue4"
col_cell_plate[targets$Plaat.nr ==  "A7"] = "deepskyblue2"
col_cell_plate[targets$Plaat.nr ==  "A8"] = "cyan"
col_cell_plate[targets$Plaat.nr ==  "B1"] = "darkgreen"
col_cell_plate[targets$Plaat.nr ==  "B2"] = "forestgreen"
col_cell_plate[targets$Plaat.nr ==  "B3"] = "chartreuse4"
col_cell_plate[targets$Plaat.nr ==  "B4"] = "chartreuse3"
col_cell_plate[targets$Plaat.nr ==  "B5"] = "chartreuse2"
col_cell_plate[targets$Plaat.nr ==  "B6"] = "chartreuse"
col_cell_plate[targets$Plaat.nr ==  "B7"] = "green3"
col_cell_plate[targets$Plaat.nr ==  "B8"] = "green1"
col_cell_plate[targets$Plaat.nr ==  "C1"] = "firebrick4"
col_cell_plate[targets$Plaat.nr ==  "C2"] = "firebrick2"
col_cell_plate[targets$Plaat.nr ==  "C3"] = "brown2"
col_cell_plate[targets$Plaat.nr ==  "C4"] = "brown1"
col_cell_plate[targets$Plaat.nr ==  "C5"] = "tomato"
col_cell_plate[targets$Plaat.nr ==  "C6"] = "brown2"
col_cell_plate[targets$Plaat.nr ==  "C7"] = "orangered2"
col_cell_plate[targets$Plaat.nr ==  "C8"] = "red"
pdf("../Plots_APP23_RNASEQ/MDS_plot_plate_number.pdf") 
plotMDS.DGEList(dge, col= col_cell_plate, main="MDS plot (plot number effect)")
dev.off() 

####################################################################
#           Plot for the plate (focussing on the number)           #
####################################################################
col_cell_plate_nr <- rep("black", dim(targets)[1])
col_cell_plate_nr[targets$Plaat.nr ==  "A1"] = "brown4"
col_cell_plate_nr[targets$Plaat.nr ==  "A2"] = "red"
col_cell_plate_nr[targets$Plaat.nr ==  "A3"] = "orange"
col_cell_plate_nr[targets$Plaat.nr ==  "A4"] = "yellow"
col_cell_plate_nr[targets$Plaat.nr ==  "A5"] = "green"
col_cell_plate_nr[targets$Plaat.nr ==  "A6"] = "blue"
col_cell_plate_nr[targets$Plaat.nr ==  "A7"] = "magenta"
col_cell_plate_nr[targets$Plaat.nr ==  "A8"] = "grey"
col_cell_plate_nr[targets$Plaat.nr ==  "B1"] = "brown4"
col_cell_plate_nr[targets$Plaat.nr ==  "B2"] = "red"
col_cell_plate_nr[targets$Plaat.nr ==  "B3"] = "orange"
col_cell_plate_nr[targets$Plaat.nr ==  "B4"] = "yellow"
col_cell_plate_nr[targets$Plaat.nr ==  "B5"] = "green"
col_cell_plate_nr[targets$Plaat.nr ==  "B6"] = "blue"
col_cell_plate_nr[targets$Plaat.nr ==  "B7"] = "magenta"
col_cell_plate_nr[targets$Plaat.nr ==  "B8"] = "grey"
col_cell_plate_nr[targets$Plaat.nr ==  "C1"] = "brown4"
col_cell_plate_nr[targets$Plaat.nr ==  "C2"] = "red"
col_cell_plate_nr[targets$Plaat.nr ==  "C3"] = "orange"
col_cell_plate_nr[targets$Plaat.nr ==  "C4"] = "yellow"
col_cell_plate_nr[targets$Plaat.nr ==  "C5"] = "green"
col_cell_plate_nr[targets$Plaat.nr ==  "C6"] = "blue"
col_cell_plate_nr[targets$Plaat.nr ==  "C7"] = "magenta"
col_cell_plate_nr[targets$Plaat.nr ==  "C8"] = "grey"
pdf("../Plots_APP23_RNASEQ/MDS_plot_plate_number2.pdf") 
plotMDS.DGEList(dge, col= col_cell_plate_nr, main="MDS plot (plot number effect by letter)")
dev.off() 

####################################################################
#           Plot for the plate (focussing on the letter)           #
####################################################################
col_cell_plate_letter <- rep("black", dim(targets)[1])
col_cell_plate_letter[targets$Plaat.nr ==  "A1"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A2"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A3"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A4"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A5"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A6"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A7"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "A8"] = "cyan"
col_cell_plate_letter[targets$Plaat.nr ==  "B1"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B2"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B3"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B4"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B5"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B6"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B7"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "B8"] = "green"
col_cell_plate_letter[targets$Plaat.nr ==  "C1"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C2"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C3"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C4"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C5"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C6"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C7"] = "red"
col_cell_plate_letter[targets$Plaat.nr ==  "C8"] = "red"
pdf("../Plots_APP23_RNASEQ/MDS_plot_plate_number3.pdf") 
plotMDS.DGEList(dge, col= col_cell_plate_letter, main="MDS plot (plot number effect by number)")
dev.off() 

####################################################################
#                      Plot for the running lane                   #
####################################################################
col_cell_running_lane <- rep("black", dim(targets)[1])
col_cell_running_lane[targets$Run..lane. ==  "1"] = "red"
col_cell_running_lane[targets$Run..lane. ==  "2"] = "blue"
pdf("../Plots_APP23_RNASEQ/MDS_plot_running_lane.pdf") 
plotMDS.DGEList(dge, col= col_cell_running_lane, main="MDS plot (plot running lane)")
dev.off() 