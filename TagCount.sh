#!/usr/bin/env bash

# written by C. Kenkel, June 2020
# counts number of high-quality (MAPQ>=23, which
# corresponds to no more than 3 mismatches at Q20))
# 2bRAD reads mapping to each reference in a 
# supplied species reference list
# By default this script discards tags exhibiting
# equal match scores to multiple reference tags (XS:i:0)
# Note that species names must be appended to tags in
# the reference itself
# Note that MAPQ thresholds can also be altered by 
# making changes to the awk command within the script
# usage: TagCount.sh samfile REFlist 

SAM=$1
SppName=$2

if [ "$#" -ne 2 ]
	then
	echo ""
	echo "Usage: TagCount.sh samfile REFlist"
	echo "Where"
	echo "	samfile:	.sam file mapped against combined species reference"
	echo "	REFlist:	a list of species names included in the reference"
	echo "For example. 'TagCount.sh PN38.sam SppName'"
	echo ""
	exit
fi


for name in `cat $SppName`
do
    echo -e "$name\t" | tr -d '\n'
    val= grep -v 'XS:i:0' $SAM | awk '$5 >= 23' | grep -v '@SQ' | grep -c $name
    echo $val | tr -d '\n'
done

