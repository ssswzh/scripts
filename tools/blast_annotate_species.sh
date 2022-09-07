#!/usr/bin/sh
# @Author: Siwen Zhang, zhangsiwen@grandomics.com, June 19 2019 
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
    then makeblastdb -in ${db} -dbtype nucl -title ${db} -out ${db}
    fi
else echo "Database not found"
    exit 
fi

printf "blastn -db ${db} -query ${fasta} -evalue 1e-5 -outfmt 7 -num_threads 8 -out ${outfile}.blast \n\n"

printf "grep -v \"#\" ${outfile}.blast|awk 'BEGIN{FS=\"\\\t\";OFS=\"\\\t\"}{if(a[\$1]){if(a[\$1]<\$12){a[\$1]=\$12;b[\$1]=\$2;c[\$1]=\$0}}else{a[\$1]=\$12;b[\$1]=\$2;c[\$1]=\$0}}END{for(i in a){print c[i]}}' > ${outfile}.blast.besthit \n\n"

#printf "awk 'NR==FNR{a[\$1]=\$0} NR>FNR{if(a[\">\"\$2]){print \$0\"\\\t\"a[\"\>\"\$2]} else{print \$0\"\\\tNot found\"} }' /export/home/zhangsiwen/database/nt/nt_20170615.fa.species ${outfile}.blast.besthit > ${outfile}.blast.besthit.species \n"

#printf "cut -f13 ${outfile}.blast.besthit.species|cut -d' ' -f2-3|sed 's/ /_/g'|sed 's/[,]//g'|sort|uniq -c|sort -k1,1nr|awk '{print $2"\t"$1}' > ${outfile}.blast.besthit.species.counts \n"
