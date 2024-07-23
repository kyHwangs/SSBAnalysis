#! /bin/sh

mkdir ROOT

for FILE in ./*
do
	if [ -d $FILE ]
	then
		echo ./ROOT/${FILE}.root
		hadd ./ROOT/${FILE}.root ${FILE}/*.root	
	fi
done


rm ./ROOT/ROOT.root

hadd ./ROOT/UL2016APV_Data.root ./ROOT/UL2016APV_Run2016*.root
hadd ./ROOT/UL2016_Data.root ./ROOT/UL2016_Run2016*.root
hadd ./ROOT/UL2017_Data.root ./ROOT/UL2017_Run2017*.root
hadd ./ROOT/UL2018_Data.root ./ROOT/UL2018_Run2018*.root

