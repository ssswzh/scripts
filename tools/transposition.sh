#!/usr/bin/bash

if [ $1 == "-h" -o $1 == "--help" ]
then printf "\n"
     printf "Usage:\n"
     printf "sh $0 file_to_be_transposed (stdout) \n\n"
     exit 1
fi

infile=$1

awk '
BEGIN{FS="\t"}
{ 
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {    
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str"\t"a[i,j];
        }
        print str
    }
}' $infile 
