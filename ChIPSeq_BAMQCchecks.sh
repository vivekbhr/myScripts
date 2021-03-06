#!/bin/bash
###

echo "This file is for qc checks of the bam files for ChIP-Seq. 
It includes bamCorrelate, computeGCbias and bamFingerprint"

echo "
Input 1 : Path to folder with all the bam files to analyse, along with the bam indexs
Input 2 : bamcorrelate- correlation method
Input 3 : output folder path (new directories will be created here)
Input 4 : scriptname (bamC, GCbias, bamF, all)
"
## This script should itself create directories and place the files there

#####--------------
#export PYTHONPATH=/data/projects/akhtar/bhardwaj/project1/Gitrepo/deepTools:$PYTHONPATH
scriptfolder=/package/deepTools-1.5.9.1/bin/
workfolder=${3}
infolder=${1}
scriptname=${4}

#-------------------

if [ ${4} == "bamC" ]
then
	## Correlate between different female files
	echo "bamcorrelate Working"

	mkdir ${workfolder}/bamcorrelate_output
	${scriptfolder}/bamCorrelate bins --corMethod ${2} --colorMap Greens \
	-b ${infolder}/*.bam -o ${workfolder}/bamcorrelate_output/bamCorr_${2}_bamcorrelate.png &


elif [ ${4} == "GCbias" ]
then

	## Compute GC bias
	echo "computing GC bias"

	mkdir ${workfolder}/gcbias_output
	for file in ${infolder}/*.bam
	do
	name=$(basename "$file")
	${scriptfolder}/computeGCBias --bamfile ${file} --effectiveGenomeSize 2150570000 --genome /data/projects/misc/genomes/2bit/mm9.2bit \
	--fragmentLength 300 -freq ${workfolder}/gcbias_output/${name}_GCfreq_file --biasPlot ${workfolder}/gcbias_output/${name}_GCplot.png &
	done

elif [ ${4} == "bamF" ]
then
	## bamFingerprint

	echo "bamFingerprint Working"
	mkdir ${workfolder}/bamFP_output
	for file in ${infolder}/*.bam
	do
	name=$(basename "$file" | sed 's/.//21g')
	${scriptfolder}/bamFingerprint --bamfiles ${infolder}/InputA_femES_129S1.CASTEiJ_suspendersMerged.bam \
	${infolder}/InputB_femES_129S1.CASTEiJ_suspendersMerged.bam \
	${infolder}/InputA_femNPC_129S1.CASTEiJ_suspendersMerged.bam \
	${infolder}/InputB_femNPC_129S1.CASTEiJ_suspendersMerged.bam \
	${file} \
	--plotFile ${workfolder}/bamFP_output/${name}_bamFP_plot.png --labels  --outRawCounts ${workfolder}/bamFP_output/${name}_outRawCounts &
	done


elif [ ${4} == "all" ] ## Do all

then
	echo "bamcorrelate Working"

	mkdir ${workfolder}/bamcorrelate_output
	${scriptfolder}/bamCorrelate bins --corMethod ${2} --colorMap Greens \
	-b ${infolder}/*.bam -o ${workfolder}/bamcorrelate_output/bamCorr_${2}_bamcorrelate.png &


	#-----------
	echo "computing GC bias"

	mkdir ${workfolder}/gcbias_output
	for file in ${infolder}/*.bam
	do
	name=$(basename "$file")
	${scriptfolder}/computeGCBias --bamfile ${file} --effectiveGenomeSize 2150570000 --genome /data/projects/misc/genomes/2bit/mm9.2bit \
	--fragmentLength 300 -freq ${workfolder}/gcbias_output/${name}_GCfreq_file --biasPlot ${workfolder}/gcbias_output/${name}_GCplot.png &
	done
	#-----------
	echo "bamFingerprint Working"

	mkdir ${workfolder}/bamFP_output
	for file in ${infolder}/*.bam
	do	
	name=$(basename "$file" | sed 's/.//21g')
	${scriptfolder}/bamFingerprint --bamfiles ${infolder}/$(ls ${infolder} | grep dupremoved.bam$ | grep InputA | head -1 | tail -1) \
	${infolder}/$(ls ${infolder} | grep dupremoved.bam$ | grep InputB | head -1 | tail -1) \
	${infolder}/$(ls ${infolder} | grep dupremoved.bam$ | grep InputA | head -2 | tail -1) \
	${infolder}/$(ls ${infolder} | grep dupremoved.bam$ | grep InputB | head -2 | tail -1) \
	${file} \
	--plotFile ${workfolder}/bamFP_output/${name}_bamFP_plot.png --labels InputA_femES InputB_femES InputA_femNPC InputB_femNPC ${name} --outRawCounts ${workfolder}/bamFP_output/${name}_outRawCounts &
	done

#---------------
else echo "Please give a scriptname"
fi
