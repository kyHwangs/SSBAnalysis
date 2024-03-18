#! /bin/sh

PROCESS=$1
TEMP=`ll ../../input/$PROCESS/ | grep list | wc -l`
ANUM=1
LISTNUM=$(($TEMP - $ANUM))
LOOP=$(seq 1 $LISTNUM)

MISSING_ROOT=()

for i in $LOOP
do
	if [ -f ../../input/$PROCESS/$PROCESS\_$i.list ];
	then
		if [ -f ./$PROCESS/$PROCESS\_$i.root ];
		then
			echo $i LIST ROOT EXIST
		else
			echo $i LIST      EXIST
			MISSING_ROOT+=($i)
		fi
	else
		if [ -f ./$PROCESS/$PROCESS\_$i.root ];
		then
			echo $i ROOT      EXIST
		else
			echo $i NO BOTH FILE
			MISSING_ROOT+=($i)
		fi
	fi
done

echo MISSING ROOT Files : $MISSING_ROOT

