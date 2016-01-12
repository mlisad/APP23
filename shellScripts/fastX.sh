#!/bin/bash
#
# The first step is to download the fastx toolkit.
# http://hannonlab.cshl.edu/fastx_toolkit/download.html
#
# This file is used to perform the trimming on the fastQ files.
# the length of the trimming can be determined fy the FastQC.
# The input that must be given is the length of trimming (example: 68)

############################################################
#					      Pathways						   #
############################################################
# The pathways can be changed to the pathway of your own.
fastXPackagePathway="/home/mdubbelaar/Desktop/AligningAndQC/FastX/bin"
pathway="/media/mdubbelaar/BD_4T"
# The direcory created can be used to store the outcome of the fastX.
mkdir -p $pathway/MouseADSamples/
############################################################
#							Code						   #
############################################################
# For each fastQ file in the fastQFilesPathway
# Trim the dataset, save it in the outputfile, trim it with the length 
# that is profided into the commandline and run it in the background. 
for fastq in $pathway/FastQ/*.fastq;
do
	# Takes the 6th element of a spliced fastq string, which is the filename.
	sample="$(echo $fastq | awk -F '[/]' '{print $6}')"
	# Beware the -Q33 is added since the data is analysed with the use of illumina.
	$fastXPackagePathway/fastx_trimmer -i $fastq  -o $pathway/MouseADSamples/$sample -l $1 -Q33 &
# Forloop is closed
done
