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
hadd ./ROOT/Data_Run2016_singleMu.root ./ROOT/Data_SingleMuon_*.root

