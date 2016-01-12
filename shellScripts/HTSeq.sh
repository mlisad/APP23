#!/bin/bash
#
# This script will create a .count data with the results that are 
# obtained from the files that are created before running the RNA SeQC.
#
# Make sure that you install HTSeq correctly.
#
# Install HTSeq by terminal 
# (sudo apt-get install python-htseq)
# Install other files needed for HTSeq 
# (sudo apt-get install build-essential python2.7-dev python-numpy python-matplotlib)
<<<<<<< HEAD
############################################################
#					      Pathways						   #
############################################################
# Makes sure that the string which is given when calling the file 
# is converted into a capitalized string
=======
#

# The right gtfFile is saved as gtfFile according to the organism.
# Also the directory is chosen this way.
>>>>>>> bbb7012ddba2c01276e0db230fcb806bee3a8e85
organism="$(echo $1 | tr [A-Z] [a-z] | sed -e 's/^./\U&/g; s/ ./\U&/g')"
pathwayStart=/media/mdubbelaar/BD_4T
# If the organism == Mouse or Human it will be used to make sure that
# the gtfFile is saved correctly. 
if [ "$organism" == "Mouse" ]; then
	gtfFile=$pathwayStart/Mus_musculus.NCBIM37.67.gtf
elif [ "$organism" == "Human" ]; then
	gtfFile=$pathwayStart/Homo_sapiens.GRCh37.71.gtf
else
	printf "\033[0;31mPlease define mouse or human when executing this script\n\033[0m"
	exit
# finish the if-else statement
fi
<<<<<<< HEAD
############################################################
#							Code						   #
############################################################
for bamFile in $pathwayStart/APP23/RNA-SeqC/$organism/SamSORTED-*.bam; 
do 	
	# Filters the name from the file and uses it for the count file.
	file="$(echo $bamFile | awk -F '[/]' '{print $8}')"
	htseq-count -f bam $bamFile -s no -r name $gtfFile > '/media/mdubbelaar/BD_4T/APP23/HTSeq/'$organism'/'$file'.count'
# Forloop is closed
=======
# The counts of each gene are saved into a .counts text file.
# This process is done for each samsorted bam file that is saved into the directory. 
for i in SamSORTED-*.bam; 
do 
	htseq-count -f bam $i -s no -r name $gtfFile > '/media/mdubbelaar/BD_4T/APP23/HTSeq/'$organism'/'$i'.count'
>>>>>>> bbb7012ddba2c01276e0db230fcb806bee3a8e85
done

