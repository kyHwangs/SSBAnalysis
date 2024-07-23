#! /bin/sh

for FILE in ./*
do
	if [ -d $FILE ]
	then
		ROOTNUM=`ll $FILE | grep root | wc -l`
		LISTNUM=`ll ../../input/$FILE/ | grep list | wc -l`
		echo $FILE $ROOTNUM $LISTNUM $(($LISTNUM-$ROOTNUM))
	fi
done
