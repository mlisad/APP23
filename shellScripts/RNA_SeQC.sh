#!/bin/bash
# Execute the file with bash ./RNA_SeQC.sh the '[]' are not includes in the sh functionality
# Mouse of Human must be defined when this script is called.
#
# The data from STAR will be sorted with picard
# The readgroups will be replaces or added
# The groups will be sorted
# The bam file will be reorderd
# The bam file will be indexed 
# A RNA-SeQC will be performed to determine the quality check of the data.
#
# Make sure that a sampleFile.csv is created as followed
# sampleID\t"PathwayToFile/SamSORTED-$sample-GROUPS.bam"\tDescription
# This file had a header as well.
#
# Make sure to download Samtools: http://www.htslib.org/download/
# Make sure to download RNA SeQC: https://www.broadinstitute.org/cancer/cga/rnaseqc_download
# Make sure to download Picard:
# wget https://github.com/broadinstitute/picard/releases/download/1.140/picard-tools-1.140.zip -O picard-tools-1.140.zip
# unzip picard-tools-1.140.zip

############################################################
#					      Pathways						   #
############################################################
pathwayStart=/media/mdubbelaar/BD_4T
picard=/home/mdubbelaar/Desktop/AligningAndQC/afterStar/Picard/picard-tools-1.140/picard.jar
RNASeQC=/home/mdubbelaar/Desktop/AligningAndQC/afterStar/RNA_SeqC/RNA-SeQC_v1.1.8.jar
# Makes sure that the string which is given when calling the file 
# is converted into a capitalized string.
organism="$(echo $1 | tr [A-Z] [a-z] | sed -e 's/^./\U&/g; s/ ./\U&/g')"

# Creates the directory where the RNA-SeqC data is stored
mkdir -p $pathwayStart/APP23/RNA-SeqC/$organism
rnaSeqcPathway=$pathwayStart/APP23/RNA-SeqC/$organism

# If the organism == Mouse or Human it will be used to create the
# fasta.fai and the .dict file of the fasta. The correct fastaFile 
# pathway, gtfFile and sampleFile are given here as well.
if [ "$organism" == "Mouse" ]; then
	fastaFile=$pathwayStart/Mus_Musculus.NCBIM37.67.fasta
	gtfFile=$pathwayStart/Mus_musculus.NCBIM37.67.gtf
	sampleFile=$pathwayStart/MiceSampleFile.csv
	# The underlying code must not be forgotten!!
	# The reference which is used for picard needs a dict along with it.
	printf "\033[0;36mCreating the SequenceDictionary and the fasta reference\n\033[0m"
	echo -e java -jar $picard CreateSequenceDictionary REFERENCE=$fastaFile OUTPUT=$pathwayStart/Mus_Musculus.NCBIM37.67.dict
	printf "\033[0;36mCreating the fasta reference\n\033[0m"
	echo -e samtools faidx $fastaFile
	# ------------------------------------------------------------------
elif [ "$organism" == "Human" ]; then
	gtfFile=$pathwayStart/Homo_sapiens.GRCh37.71.gtf
	fastaFile=$pathwayStart/human_g1k_v37.fasta
	sampleFile=$pathwayStart/HumanSampleFile.csv
	# The underlying code must not be forgotten!!
	# The reference which is used for picard needs a dict along with it.
	printf "\033[0;36mCreating the SequenceDictionary\n\033[0m"
	java -jar $picard CreateSequenceDictionary R=$fastaFile O=$pathwayStart/human_g1k_v37.dict
	printf "\033[0;36mCreating the fasta reference\n\033[0m"
	samtools faidx $fastaFile
	# ------------------------------------------------------------------
else
	# If the organism is not equal to Mouse or Human an error will be given
	# an the file will exit.
	printf "\033[0;31mPlease define mouse or human when executing this script\n\033[0m"
	exit
# finish the if-else statement
fi
############################################################
#							Code						   #
############################################################
for i in $pathwayStart/STAR_Align/$organism/*/*.sam;
do 
	# Getting the filename to make sure that the output file gets an unique name
	sample="$(echo $i | awk -F '[/]' '{print $7}')"
	# These directories below are made to make sure that the files that are created for the samtools sortsset.
	# The files are temporarely stored into the tmp file which will be deleted later.
	mkdir -p $rnaSeqcPathway/tmp/
	mkdir -p $rnaSeqcPathway/tmp/$sample
	# The first step is to sort the retreived sam file, this file is made with the use of the samsort. 
	printf "\033[0;36mBusy for sample $sample \n\033[0m"
	printf "\033[0;32mSorting\n\033[0m"
	java -jar $picard SortSam I=./$i O=$rnaSeqcPathway/$sample.bam SO=queryname
	# The second step to add or replace read groups.
	printf "\033[0;32mAdding or replacing the read groups\n\033[0m"
	java -jar $picard AddOrReplaceReadGroups I=$rnaSeqcPathway/$sample.bam O=$rnaSeqcPathway/$sample-GROUPS.bam LB=whatever PL=illumina PU=whatever SM=whatever
	# The third step reorders the bam file.
	printf "\033[0;32mReordering\n\033[0m"
	java -jar $picard ReorderSam I=$rnaSeqcPathway/$sample-GROUPS.bam O=$rnaSeqcPathway/SORT-$sample-GROUPS.bam R=$fastaFile
	# Trying to sort the dataset again cause the data needs to be sorted before indexing.
	printf "\033[0;32mSam sorting\n\033[0m"
	samtools sort -f $rnaSeqcPathway/SORT-$sample-GROUPS.bam -o $rnaSeqcPathway/tmp/$sample/SamSORTED-$sample-GROUPS.bam
	# Because of some difficulties, the samtools sort, doesnt merge by itself, that is why the merge function is called by its own.
	printf "\033[0;32mSam Merging\n\033[0m"
	samtools merge $rnaSeqcPathway/SamSORTED-$sample-GROUPS.bam $rnaSeqcPathway/tmp/$sample/SamSORTED-*.bam
	# After the reordering the files must be sorted, the first sorting is probably not good enough to sort the items in the bam file properly.
	printf "\033[0;32mIndexing\n\033[0m"
	samtools index $rnaSeqcPathway/SamSORTED-$sample-GROUPS.bam
	# The unnecessary files and directories are deleted again.
	printf "\033[0;32mRemoving unnecessary files\n\033[0m"
	rm $rnaSeqcPathway/tmp/$sample/*
	rmdir $rnaSeqcPathway/tmp/$sample/
	rmdir $rnaSeqcPathway/tmp/
	rm $rnaSeqcPathway/SORT*
	rm $rnaSeqcPathway/$sample-GROUPS.bam
	done
# RNA_SeqC (This steps takes a couple of hours)
# This step makes sure that the quality of the alignment is measured.
# The results of this Quality Check can be found in the RNA-SeQC_Report
printf "\033[0;32mPerforming the RNA_SeQC on sample $sample\n\033[0m"
java -jar $RNASeQC -s $sampleFile -t $gtfFile -r $fastaFile -o rnaSeqcPathway/RNA-SeQC'_'Report
