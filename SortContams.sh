#!/usr/bin/env bash

# written by C. Kenkel, July 2020
# Identifies contaminant sequences post clustering
# (using cd-hit-est) based on taxonomic assignment 
# resulting from blastn match of
# all original sequences to taxonomic database
# e.g. NCBI nt or similar. 
# Note that cd-hit-est may truncate sequence IDs in
# the .clstr output file. If so, the Contams.list 
# will need to be correspondingly truncated for this 
# script to work
# the script takes the number of the sequence from the .clstr
# file and subtracts it from the line number of the file
# to return the reference sequence from cd-hit-est
# this is needed because while the reference may not be 
# identified as a contaminant, a closely related variant
# could be
# usage: SortContams.sh Contams.list Ref.fasta.clstr Ref.fasta

LIST=$1
CLSTR=$2
FASTA=$3

if [ "$#" -ne 3 ]
        then
	echo ""
        echo "Usage: SortContams.sh Contams.list Ref.fasta.clstr Ref.fasta"
        echo "Where"
        echo "  Contams.list:   a list of known contaminants, truncated if necessary to match .clstr output"
        echo "  Ref.fasta.clstr:the .clstr output file from cd-hit-est identifying similar sequence clusters"
        echo "  Ref.fasta:      the reduced .fasta output file from the cd-hit-est clustering"
        echo "For example. 'SortContams.sh Contams2remove_trunc.tab Paus_bcgIN1sam_clstr.fasta.clstr Paus_bcgIN1sam_clstr.fasta'"
        echo ""
        exit
fi

exec < $LIST 
while read LINE 
do
    seqnum="$(grep "$LINE" $CLSTR | cut -f1)"
    lineNum="$(grep -n "$LINE" $CLSTR | head -n 1 | cut -d: -f1)"
    NumToPrint="$(($lineNum - $seqnum))"
    seq="$(awk "NR==$NumToPrint{print;exit}" $CLSTR | cut -d' ' -f 2 | cut -d'_' -f 1,2,3 | sed 's/>//g')"
    grep "$seq" $FASTA | cut -d' ' -f 1 | sed 's/>//g'
#    echo $seqnum
#    echo $lineNum
#    echo $NumToPrint
#    echo $seq
done

