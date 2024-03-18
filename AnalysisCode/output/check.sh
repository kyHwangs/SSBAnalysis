#! /bin/sh

for FILE in ./*
do
	if [ -d $FILE ]
	then
		echo $FILE
		ROOTNUM=`ll $FILE | grep root | wc -l`
		LISTNUM=`ll ../../input/$FILE/ | grep list | wc -l`
		echo $ROOTNUM $LISTNUM $(($ROOTNUM-$LISTNUM))
	fi
done
