#!/usr/bin/bash

# Created by Zhang Siwen
# History:
#     20200115, manuscript


if [ $# -lt 2 -o $# -gt 3 ];
    then
    printf "\nUsage: sh $0 <sample_list> <outdir> <optional: palindrome_flank_size> \n\n";
    printf "Formats:\n"
    printf "sample_list: sampleprefix, chr, start, end, vcf, bam [[File, Tab-delimited, without header]] \n"
    printf "        e.g. sample1    X   31882243    32077759    /path/to/vcf    /path/to/bam \n"
    printf "        default ref: /export/database/pub_database/hg19/human_g1k_v37.fasta \n"
    printf "outdir: out path, both absolute path and relative path \n"
    printf "palindrome_flank_size: flank size near breakpoint to find palindrome \n\n"
    exit 1;
fi

# Arguments
samplelist=$1
outdir=$2
flanksize=$3
ref="/export/home/zhangsiwen/database/hg19/hg19.fa"


if [ ! -d ${outdir} ]
then mkdir -p ${outdir}
fi


while read sample chr start end vcf bam
do

# get readid from vcf
awk -v chrom=${chr} -v start=${start} '{if($1==chrom&&$2==start){print $8}}' ${vcf}|cut -d';' -f10|sed 's/RNAMES=//g'|sed 's/,/\n/g' > ${outdir}/${sample}.${chr}_${start}_${end}.readid

# extract region
samtools view -h ${bam} > ${outdir}/${sample}.sam
awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR{a[$1]=$1}NR>FNR{if($1~/^@/){print $0} if(a[$1]){print $0}}' ${outdir}/${sample}.${chr}_${start}_${end}.readid ${outdir}/${sample}.sam > ${outdir}/${sample}.${chr}_${start}_${end}.readid.sam
rm ${outdir}/${sample}.sam
#cat <(samtools view -h ${bam}) ${outdir}/${sample}.${chr}_${start}_${end}.readid.aln > ${outdir}/${sample}.${chr}_${start}_${end}.readid.sam
samtools view -Sb ${outdir}/${sample}.${chr}_${start}_${end}.readid.sam > ${outdir}/${sample}.${chr}_${start}_${end}.bam
samtools index ${outdir}/${sample}.${chr}_${start}_${end}.bam

# decode cigar
python /export/home/zhangsiwen/scripts/tools/print_aln.py -r ${ref} -i ${outdir}/${sample}.${chr}_${start}_${end}.bam > ${outdir}/${sample}.${chr}_${start}_${end}.bam.aln

# decode mapping position
python /export/home/zhangsiwen/scripts/tools/decode_mapping_position.py --bam ${outdir}/${sample}.${chr}_${start}_${end}.bam --out ${outdir}/${sample}.${chr}_${start}_${end}.pos --seq True --plot False

# get base between splits
awk 'BEGIN{FS="\t";OFS="\t"}NR==FNR{a[$1]=$1}NR>FNR{n=0;if(a[$1]){if(n<3){print $0; n++}else{a[$1]=0}}}' ${outdir}/${sample}.${chr}_${start}_${end}.readid ${outdir}/${sample}.${chr}_${start}_${end}.pos|sort|paste - -|awk '{print $1"\t"$2"\t"$3"\t"$13"\t"$14"\t"$22}'|awk '{split($2","$3","$4","$5,a,","); asort(a); printf($1"\t"); for(i=1;i<=length(a);i++)printf("%s\t",a[i]); printf($6"\n");}'|awk '{base=substr($6,$3+1,$4-$3); print $0"\t"base}' > ${outdir}/${sample}.${chr}_${start}_${end}.bam.pos.sort_pos
cut -f7 ${outdir}/${sample}.${chr}_${start}_${end}.bam.pos.sort_pos|sort|uniq -c|sort -k1,1nr|awk '{print $2"\t"$1}' > ${outdir}/${sample}.${chr}_${start}_${end}.bam.pos.sort_pos.count

done < ${samplelist}



if [ ${flanksize} == "" ]
    then echo "palindrome_flank_size is not set, pass the step to find palindrome around breakpoints.";
        echo "Finished."
    else
        awk -v size=${flanksize} 'BEGIN{FS="\t";OFS="\t"}{print $2,$3-size,$3+size,$1"\n"$2,$4-size,$4+size,$1}' ${samplelist} > ${outdir}/${samplelist}.breakpoint_flank${flanksize}.bed
        bedtools nuc -fi ${ref} -bed ${outdir}/${samplelist}.breakpoint_flank${flanksize}.bed -seq > ${outdir}/${samplelist}.breakpoint_flank${flanksize}.seq
        python2 /export/home/zhangsiwen/scripts/tools/find_all_palindrome.py ${outdir}/${samplelist}.breakpoint_flank${flanksize}.seq 13 
        echo "Finished."
fi
