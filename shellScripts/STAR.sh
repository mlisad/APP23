#!/bin/bash
# The file needs to be calles with (Human or Mice).
#
# First step is to download the STAR aligner https://github.com/alexdobin/STAR/releases
# 
# Depending on the organism the fasta and gtf files are use to create a genome directory
# The alignment of the dataset against that newly created genome dir will be the next step.
# The last step is optional when the other file is also known, 
# this file will perform the RNA-SeQC and the steps between the STAR alignment
# and the RNA-SeQC.

############################################################
#					      Pathways						   #
############################################################
# The pathway to the STAR programm is defined
# And the file need an input (Mouse/Human)
STAR=/home/mdubbelaar/Desktop/AligningAndQC/STAR/STAR-STAR_2.4.2a/bin/Linux_x86_64_static/STAR
pathwayStart=/media/mdubbelaar/BD_4T
# Makes sure that the string which is given when calling the file 
# is converted into a capitalized string.
organism="$(echo $1 | tr [A-Z] [a-z] | sed -e 's/^./\U&/g; s/ ./\U&/g')"

# An if-else statement will be walked through, this statement makes sure that the right files will be called
# for the human or mouse datasets. A print will be returned with the message that mouse or human must be
# defined when this script is used.
if [ "$organism" == "Mouse" ]; then
fastQPathway=$pathwayStart/MouseADSamples
fastaFile=$pathwayStart/Mus_Musculus.NCBIM37.67.fasta
gtfFile=$pathwayStart/Mus_musculus.NCBIM37.67.gtf
elif [ "$organism" == "Human" ]; then
fastQPathway=$pathwayStart/HumanADSamples
fastaFile=$pathwayStart/human_g1k_v37.fasta
gtfFile=$pathwayStart/Homo_sapiens.GRCh37.71.gtf
else
	printf "\033[0;31mPlease define mouse or human when executing this script\n\033[0m"
	exit
# finish the if-else statement
fi
############################################################
#							Code						   #
############################################################
# Creates a directory where the genome can be generated.
mkdir -p $pathwayStart/StarGenomeGenerated/$organism
# The first step is to create the genome with the use of the fasta file and the gtf file of the right species.
# After the creation of the genome, the dataset will be aligned against this genome.
$STAR --runMode genomeGenerate --genomeDir $pathwayStart/StarGenomeGenerated/$organism --genomeFastaFiles $fastaFile --sjdbGTFfile $gtfFile --runThreadN 5
# The items in the given directory will be walked through
for fastq in $fastQPathway/*.fastq;
do
	# For each element within the file the following sentences will be returned
	echo -e "Started Alignment on $fastq\n"
	# Takes the 6th element of a spliced fastq string, which is the filename.
	file="$(echo $fastq | awk -F '[/]' '{print $6}')"
	# A directory is created to find the files in their own folder and the Star alignment is done for each file.
	mkdir -p $pathwayStart/STAR_Align/$organism/${file/.*/}/
	$STAR --genomeDir $pathwayStart/StarGenomeGenerated/$organism --runThreadN 4 --readFilesIn $fastq --outFileNamePrefix $pathwayStart/STAR_Align/$organism/${file/.*/}/ --outFilterMultimapNmax 1 --outFilterMismatchNmax 2
# The for-loop is closed
done
############################################################
#					  	  Continue						   #
############################################################
# It is optional to call another script in this script.
# The script below performs the RNA-SeQC and the steps that needs to be done before that.
. /home/mdubbelaar/shellScripts/RNA_SeQC.sh $organism
