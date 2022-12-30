#!/usr/bin/sh
# @Author: Siwen Zhang, zhangsiwen, June 19 2019 
# 20220726, check DB, add DB

if [ $# -lt 1 -o $# -ne 3 ];
    then
    echo "Usage: sh $0 <db> <raw fasta> <outprefix> > run.sh";
    exit 1;
fi

# Arguments
db=$1
fasta=$2
outfile=$3

if [ -f ${db} ]
then dir=`dirname ${db}`
    if [ ! -f ${db}.nin ]
    then echo "makeblastdb -in ${db} -dbtype nucl -title ${db} -out ${db}"
    fi
else if [ ! -f ${BLASTDB}/${db} ]
        then echo "Database not found in both `pwd` and BLASTDB=${BLASTDB}"
             echo "makeblastdb -in ${db} -dbtype nucl -title ${db} -out ${db}"
             exit 
    fi
fi

printf "blastn -db ${db} -query ${fasta} -evalue 1e-5 -outfmt 7 -num_threads 8 -out ${outfile}.blast \n\n"

printf "grep -v \"#\" ${outfile}.blast|awk 'BEGIN{FS=\"\\\t\";OFS=\"\\\t\"}{if(a[\$1]){if(a[\$1]<\$12){a[\$1]=\$12;b[\$1]=\$2;c[\$1]=\$0}}else{a[\$1]=\$12;b[\$1]=\$2;c[\$1]=\$0}}END{for(i in a){print c[i]}}' > ${outfile}.blast.besthit \n\n"

if [ -f ${db}.species ]
then
    printf "awk 'NR==FNR{if(\$1~/^>/){gsub(\">\",\"\",\$1)} a[\$1]=\$0} NR>FNR{if(a[\$2]){print \$0\"\\\t\"a[\$2]} else{print \$0\"\\\tNot found\"} }' ${db}.species ${outfile}.blast.besthit > ${outfile}.blast.besthit.species \n\n"
    printf "cut -f13- ${outfile}.blast.besthit.species|sort|uniq -c|sort -k1,1nr|awk '{printf \$2\"\\\t\";for(i=3;i<=NF;i++){if(i==NF){print \$NF\"\\\t\"\$1}else{printf \$i\" \"}}}' > ${outfile}.blast.besthit.species.counts \n"
fi
if [ -f ${BLASTDB}/${db} ]
then
    printf "awk 'NR==FNR{if(\$1~/^>/){gsub(\">\",\"\",\$1)} a[\$1]=\$0} NR>FNR{if(a["\$2]){print \$0\"\\\t\"a[\$2]} else{print \$0\"\\\tNot found\"} }' ${BLASTDB}/${db}.species ${outfile}.blast.besthit > ${outfile}.blast.besthit.species \n\n"
    printf "cut -f13- ${outfile}.blast.besthit.species|sort|uniq -c|sort -k1,1nr|awk '{printf \$2\"\\\t\";for(i=3;i<=NF;i++){if(i==NF){print \$NF\"\\\t\"\$1}else{printf\$i\" \"}}}' > ${outfile}.blast.besthit.species.counts \n"
fi
