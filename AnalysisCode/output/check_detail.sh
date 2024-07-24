#! /bin/sh

for FILE in ./*
do
	if [ -d $FILE ]
	then
		# ROOTNUM=`ll $FILE | grep root | wc -l`
		LISTNUM=`ll ../../input/$FILE/ | grep list | wc -l`

		for i in `eval echo {1..$LISTNUM}`
		do
			ROOTNAME="${FILE}/${FILE}_${i}.root"
			if [ ! -f $ROOTNAME ];
			then
				echo "$ROOTNAME does not exist"
			fi
		done
	fi
done
