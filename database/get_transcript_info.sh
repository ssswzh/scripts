#!/usr/bin/bash

if [ $1 == "-h" -o $1 == "--help" -o $1 == "-help" ] || [ $# == 0 ]
then printf "\nUSAGE: \n    sh $0 hg19.ncbiRefSeq.gtf transcript.list\n\n"
printf "transcript.list format: \n"
printf "NM_004006\n"
printf "NM_004009\n\n"
exit 1
fi

gff=$1
transcriptlist=$2

for transcript in `cat $transcriptlist`
do
grep ${transcript} ${gff}|
awk '
    BEGIN{
    FS="\t";OFS="\t";t=0;n=0;
    }
    {if($3=="exon"){
    t=t+$5-$4+1; n++}
    split($9,info,";");
    split(info[2],trans," ");
    gsub("\"","",trans[2]);
    }
    END{print trans[2],t,n;}'
done 

