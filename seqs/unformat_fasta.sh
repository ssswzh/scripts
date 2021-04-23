#!/bin/bash
if [ "$1" == "" ]
then
echo "Change multi-line fasta record to one line"
echo "Usage: sh $0 fasta_file > new_fasta"
exit 1
else
awk 'NR==1{print $0}NR!=1{if($1~/^>/){print "\n"$0}else{printf $0}}' $1
fi
