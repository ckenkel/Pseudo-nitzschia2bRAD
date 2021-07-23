#!/bin/bash
# organizes a list of sequence files into sample-specific directories
# input a list of sample names

LIST=$1

if [ "$#" -ne 1 ]
	then
	echo ""
	echo "Usage: organize.sh list"
	echo "Where"
	echo "	list:	a list of sample names (typically the first part of each file name)"
	echo ""
	exit
fi

exec < "$LIST"
while read LINE
do
        echo "$LINE"
	mkdir $LINE"_dir"
	mv $LINE $LINE"_dir"
done

