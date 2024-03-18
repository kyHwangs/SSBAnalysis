#! /bin/sh

for FILE in ./*
do
	if [ -d $FILE ]
	then
		LISTNUM=`tree -C $FILE | grep list | wc -l`
		ROOTNUM=`tree -C /d0/scratch/haeun/HD/CMSSW_8_0_26_patch1/src/SSBAnalysis/AnalysisCode/output/3_Data/$FILE | grep root | wc -l`
		ROOTNUMEE=`tree -C /d0/scratch/haeun/HD/CMSSW_8_0_26_patch1/src/SSBAnalysis/AnalysisCode/output/3_MC/$FILE | grep root | wc -l`
		echo $FILE $LISTNUM $ROOTNUM $ROOTNUMEE

	fi
done
