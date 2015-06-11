#!/bin/bash

## Written by vivek.. 5-6 Nov 2014.. updated 8th Nov 2014
## Remove ENCODE blacklisted regions from the given bam file and create new bam files.. for female chip-seq data in mouse
## modify the parameters for any other cell

infolder=${1}
outfolder=${2}

#-------------------------------\
echo "This script removes ENCODE blacklisted regions and duplicated reads from the given bam file and creates new bam files..
Input 1 = folder with bam files
Input 2 = folder to put the output
Input 3 = celltype, OR any grep-readable regular expressions to pick files from input folder
input 4 = scriptname (filter,rmdup). leave blank for both
"
declare -a filename
declare -a newfilename
	totalfiles=$(ls ${infolder} | grep .bam$ | grep ${3} | wc -l) # stores the number of files in the infolder following given pattern
	for i in $(seq 1 ${totalfiles}); do filename[i]=$(ls ${infolder} | grep .bam$ | grep ${3} | head -${i} | tail -1); done # save filenames in array(starting from 1, not 0)

if [ ${4} == "filter" ]
then
	echo "Removing blacklisted and random regions"
	for num in $(seq 1 ${totalfiles})
	do
	/package/samtools-1.1/samtools view -h ${infolder}/${filename[@]:${num}:1} | python /data/akhtar/bhardwaj/Gitrepo/MultipurposeScripts_pipelines/samFilter_v3.py --filter_out_from_BED  /data/manke/repository/misc/annotations/blacklist_ENCODE/mm9-blacklist.bed --random --lowqual --chrM > ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.sam
	/package/samtools-1.1/samtools view -Sb ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.sam -o ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.bam
	done

elif [ ${4} == "rmdup" ]
then
	echo "Removing duplicate reads, assuming sorted input bam files"
	for num in $(seq 1 ${totalfiles})
	do
	java -jar /package/picard-tools-1.99/MarkDuplicates.jar REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT INPUT=${infolder}/${filename[@]:${num}:1} METRICS_FILE=${outfolder}/${filename[@]:${num}:1}.matrix OUTPUT=${outfolder}/${filename[@]:${num}:1}_dupremoved.bam
	done

else
	echo "Doing all filtering steps together (shall generate intermediate files as well)"
	echo "Removing blacklisted and random regions" # exact copy of the script in first loop
	for num in $(seq 1 ${totalfiles})
	do
	/package/samtools-1.1/samtools view -h ${infolder}/${filename[@]:${num}:1} | python /data/akhtar/bhardwaj/Gitrepo/MultipurposeScripts_pipelines/samFilter_v3.py --filter_out_from_BED /data/manke/repository/misc/annotations/blacklist_ENCODE/mm9-blacklist.bed --random --lowqual --chrM > ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.sam
	/package/samtools-1.1/samtools view -Sb ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.sam -o ${outfolder}/${filename[@]:${num}:1}_blaclist_random_removed.bam
	done

	echo "Removing duplicated reads from blacklist removed files now" # modified as compared to 2nd loop script (only takes blacklist filtered files)
	newfiles=$(ls ${outfolder} | grep blaclist_random_removed.bam$ | grep ${3} | wc -l) #new number of files (blacklist removed)
	for i in $(seq 1 ${newfiles}); do newfilename[i]=$(ls ${outfolder} | grep blaclist_random_removed.bam$ | grep ${3} | head -${i} | tail -1); done ## new array to store the filenames (blacklist removed)

if [ ! -d ${outfolder}/filtered_bam] #check if output folder for filtred files exists
then mkdir ${outfolder}/filtered_bam
else echo "Writing on the pre-created directory filtered_bam"
	fi
	for num in $(seq 1 ${newfiles})
	do
	java -jar /package/picard-tools-1.99/MarkDuplicates.jar REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT INPUT=${outfolder}/${newfilename[@]:${num}:1} METRICS_FILE=${outfolder}/filtered_bam/${newfilename[@]:${num}:1}.matrix OUTPUT=${outfolder}/filtered_bam/${newfilename[@]:${num}:1}_dupremoved.bam
	done

fi
unset filename
unset newfilename
