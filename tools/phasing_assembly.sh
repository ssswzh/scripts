#!/usr/bin/bash

# Created by Zhang Siwen, zhangsiwen@grandomics.com
# History:
#       20191210, create first version

if [ $# -lt 1 -o $# -ne 10 ];
    then
    echo "Usage: sh $0 <vcf> <bam> <outdir> <output prefix> <bed file> <pacbio|nanopore> <hg19,hg38> <thread> <assemble error rate> <medaka_model, r941_min_fast|r941_min_high|r941_prom_fast|r941_prom_high|r10_min_high|r941_min_diploid_snp> > run.sh";
    exit 1;
fi

# Arguments

vcf=$1
bam=$2
outdir=$3
outprefix=$4
bed=$5
platform=$6
ref_version=$7
thread=$8
error_rate=$9
medaka_model=${10}

if [ ${ref_version} == "hg19" ]
        then ref="/export/database/pub_database/hg19/human_g1k_v37.fasta"; 
elif [ ${ref_version} == "hg38" ]
        then ref="/export/database/pub_database/hg38_without_fragment/hg38.fa";
else echo "Reference version error."; echo "Please choose from hg19 and hg38."; exit 1;
fi


# Alignment

printf "

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=\$current_command; current_command=\$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo \"\\\\\"\${last_command}\\\" command filed with exit code \$?.\\\"' EXIT


if [ ! -d ${outdir} ]
        then mkdir ${outdir};
fi \n"

printf "

# vcf keep only snps

vcftools --vcf ${vcf} --remove-indels --recode --recode-INFO-all --out ${outdir}/${outprefix}.only_snp

bam_prefix=`basename -s .bam ${bam}`
samtools mpileup -o ${outdir}/\${bam_prefix}.mpileup -f ${ref} --output-BP --output-QNAME ${bam}

# whatshap phase

whatshap phase --changed-genotype-list ${outdir}/${outprefix}.whatshap.genotypelist.fofn --output-read-list ${outdir}/${outprefix}.whatshap.readslist.fofn --tag PS --reference ${ref} --ignore-read-groups --distrust-genotypes -o ${outdir}/${outprefix}.whatshap.vcf.gz ${outdir}/${outprefix}.only_snp.recode.vcf ${bam}

gunzip ${outdir}/${outprefix}.whatshap.vcf.gz 
\n\n"

# bed file process

printf "
# bed file process

if [ ! -d ${outdir}/assembly_canu${error_rate}_medaka${medaka_model} ]
        then mkdir ${outdir}/assembly_canu${error_rate}_medaka${medaka_model};
fi

source /export/home/zhangsiwen/.medaka_env
\n"

while read chrom start end
do
genomeSize=$(($((end))-$((start))))
printf "
# ${chrom}_${start}_${end} process starts

printf \"${chrom}\\\t${start}\\\t${end}\" > ${outdir}/region.${chrom}_${start}_${end}.bed 
bedtools intersect -a ${outdir}/${outprefix}.whatshap.vcf -b ${outdir}/region.${chrom}_${start}_${end}.bed -header > ${outdir}/${outprefix}.whatshap.${chrom}_${start}_${end}.vcf

# separate reads according to haplotype

python2 /export/home/zhangsiwen/bin/pileup2matrix.py --vcf ${outdir}/${outprefix}.whatshap.${chrom}_${start}_${end}.vcf --pileup ${outdir}/\${bam_prefix}.mpileup --prefix ${outprefix}.${chrom}_${start}_${end} --outdir ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}

# fetch reads for two haplotype

awk 'BEGIN{FS=\"\\\t\";OFS=\"\\\n\"}NR==FNR{a[\$1]=\$1}NR>FNR{if(a[\$1]){print \">\"\$1,\$10}}'  ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.list <(samtools view ${bam}) > ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.fa

awk 'BEGIN{FS=\"\\\t\";OFS=\"\\\n\"}NR==FNR{a[\$1]=\$1}NR>FNR{if(a[\$1]){print \">\"\$1,\$10}}'  ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.list <(samtools view ${bam}) > ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.fa


# canu assembly haplotype

/export/software/canu/canu-1.9/canu/Linux-amd64/bin/canu -p ${outprefix}.${chrom}_${start}_${end}.type1 -d ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.canu minReadLength=200 minOverlapLength=50 genomeSize=${genomeSize} -${platform} ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.fa -useGrid=false correctedErrorRate=${error_rate} corOutCoverage=1 stopOnLowCoverage=0

/export/software/canu/canu-1.9/canu/Linux-amd64/bin/canu -p ${outprefix}.${chrom}_${start}_${end}.type2 -d ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.canu minReadLength=200 minOverlapLength=50 genomeSize=${genomeSize} -${platform} ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.fa -useGrid=false correctedErrorRate=${error_rate} corOutCoverage=1 stopOnLowCoverage=0


# mini_assemble

/export/home/hanyue/software/pomoxis/scripts/mini_assemble -i ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.fa -o ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.mini_assemble -p ${outprefix}.${chrom}_${start}_${end}.type1.mini_assemble -t ${thread} -r ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.canu/${outprefix}.${chrom}_${start}_${end}.type1.contigs.fasta

/export/home/hanyue/software/pomoxis/scripts/mini_assemble -i ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.fa -o ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.mini_assemble -p ${outprefix}.${chrom}_${start}_${end}.type2.mini_assemble -t ${thread} -r ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.canu/${outprefix}.${chrom}_${start}_${end}.type2.contigs.fasta


# medaka_consensus

medaka_consensus -i ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.fa -d ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.mini_assemble/${outprefix}.${chrom}_${start}_${end}.type1.mini_assemble_final.fa -o ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.medaka_consensus -m ${medaka_model}
mv ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.medaka_consensus/consensus.fasta ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type1.medaka_consensus/${outprefix}.${chrom}_${start}_${end}.type1.fasta

medaka_consensus -i ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.fa -d ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.mini_assemble/${outprefix}.${chrom}_${start}_${end}.type2.mini_assemble_final.fa -o ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.medaka_consensus -m ${medaka_model}
mv ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.medaka_consensus/consensus.fasta ${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.${chrom}_${start}_${end}.type2.medaka_consensus/${outprefix}.${chrom}_${start}_${end}.type2.fasta


# ${chrom}_${start}_${end} process ends

\n" 
done < ${bed}


# result interpretation
printf "
printf \"

#***************************************
#     Phasing and assembly all done.
#***************************************

Result files:

Note: change 'region' in filename to your bed region with 'chrom_start_end'.

Two haplotype raw reads: 
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type1.fa
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type2.fa

Canu original assembled haplotype contigs:
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type1.canu/${outprefix}.region.type1.contigs.fasta
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type2.canu/${outprefix}.region.type2.contigs.fasta

Medaka_consensus results based on canu (which is final haplotype contigs):
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type1.medaka_consensus/${outprefix}.region.type1.fasta
${outdir}/assembly_canu${error_rate}_medaka${medaka_model}/${outprefix}.region.type2.medaka_consensus/${outprefix}.region.type2.fasta


\" > ${outdir}/${outprefix}.result_interpretation
"
