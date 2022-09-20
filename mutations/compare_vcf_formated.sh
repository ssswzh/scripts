#!/usr/bin/bash

# Created by Zhang Siwen
# History:
#     20200113, manuscript


if [ $# -ne 8 ];
    then
    printf "\nUsage: sh $0 <sample_list> <outdir> <baseline outprefix> <test outprefix> <exclude or include> <region.bed> <region config> <variant type> \n\n";
    printf "Formats:\n"
    printf "sample_list: SampleName, BaselineVCF(TP), TestVCF [[File, Tab-delimited]] \n"
    printf "outdir: out path, both absolute path and relative path \n"
    printf "baseline outprefix: baseline output prefix \n"
    printf "test outprefix: test output prefix \n"
    printf "exclude or include: choose from [exclude] or [include] \n"
    printf "region.bed: interested region bed to be excluded or included [[File, Tab-delimited, with header]] \n"
    printf "region config: region to be count separately, RegionName, RegionBedFile [[File, Tab-delimited]] \n"
    printf "            e.g. /export/home/zhangsiwen/project/1.DMD/DMD_compare_region.config \n"
    printf "variant type: variant type for comapring, choose from [snp], [indel] and [all]  \n"
    exit 1;
fi

# Arguments
samplelist=$1
# sampleName\ttruePosVCF\ttestVCF\ttruePosBAM\ttestBAM
outdir=$2
tppfx=$3
testpfx=$4
inex=$5
roi=$6
regionconfig=$7
varianttype=$8



# out dir
if [ ! -d ${outdir}/${varianttype}_only ]
then mkdir -p ${outdir}/${varianttype}_only 
    while read region bed
    do
    mkdir ${outdir}/${region}_${varianttype}
    done < ${regionconfig}
else echo "Error: folder exists, ${outdir}."; exit 1;
fi


# function: count snp type for each file
function count_stats(){
for i in `cat $1`
do
awk -F'_' '{print $3">"$4}' ${i}|sort|uniq -c|awk '{print $2"\t"$1}' > ${i}.count
done
}


# decide regions to be included or excluded
if [ ${inex} == "include" ]
then bedoption="--bed"
elif [ ${inex} == "exclude" ]
then bedoption="--exclude-bed"
fi

if [ "${varianttype}" == "snp" ]
then typeoption="--remove-indels"; 
elif [ "${varianttype}" == "indel" ]
then typeoption="--keep-only-indels"; 
elif [ "${varianttype}" == "all" ]
then typeoption=""; 
fi



# keep only SNPs or indels or all variants

while read sample tpvcf testvcf
do
#out=`ls ${ont}|cut -d'/' -f5|cut -d'.' -f1-3`
vcftools --vcf ${tpvcf} ${bedoption} ${roi} --out ${outdir}/${varianttype}_only/${sample}.tp_${varianttype} --recode --keep-INFO-all ${typeoption} 
grep "#\|1/1\|1|1" ${outdir}/${varianttype}_only/${sample}.tp_${varianttype}.recode.vcf|grep -v "SVTYPE=" > ${outdir}/${varianttype}_only/${sample}.${tppfx}_${varianttype}.vcf
vcftools --vcf ${testvcf} ${bedoption} ${roi} --out ${outdir}/${varianttype}_only/${sample}.test_${varianttype} --recode --keep-INFO-all ${typeoption}
grep "#\|1/1\|1|1" ${outdir}/${varianttype}_only/${sample}.test_${varianttype}.recode.vcf|grep -v "SVTYPE=" > ${outdir}/${varianttype}_only/${sample}.${testpfx}_${varianttype}.vcf
printf ${sample}"\t"${outdir}/${varianttype}_only/${sample}.${tppfx}_${varianttype}.vcf"\t"${outdir}/${varianttype}_only/${sample}.${testpfx}_${varianttype}.vcf"\n" >> ${outdir}/${varianttype}_only/vcf_file_list
done < ${samplelist}



# format_vcf.py比较差异
while read region bed
do
while read sample tpvcf testvcf
do
vcftools --vcf ${tpvcf} --bed ${bed} --out ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${varianttype}.${region} --recode --keep-INFO-all 
vcftools --vcf ${testvcf} --bed ${bed} --out ${outdir}/${region}_${varianttype}/${sample}.${testpfx}_${varianttype}.${region} --recode --keep-INFO-all 

python3 /export/home/zhangsiwen/scripts/tools/format_vcf.py --base ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${varianttype}.${region}.recode.vcf --test ${outdir}/${region}_${varianttype}/${sample}.${testpfx}_${varianttype}.${region}.recode.vcf --prefix ${sample}.${tppfx}_${testpfx}.${region}_${varianttype} --outdir ${outdir}/${region}_${varianttype}/

printf ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.tp"\t"${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.fp"\t"${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.fn"\n" >> ${outdir}/${region}_${varianttype}/vcf_file_list
done < ${outdir}/${varianttype}_only/vcf_file_list
count_stats ${outdir}/${region}_${varianttype}/vcf_file_list
done < ${regionconfig}




# 统计结果
while read region bed
do
printf "Sample\tTP\tFP\tFN\tPrecision\tRecall\tFscore\n" >> ${outdir}/${region}.${tppfx}_${testpfx}.${varianttype}.stats
while read sample tpvcf testvcf
do
result=`grep -v "Sample" ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.stats`
if [[ ${result} == "" ]]; 
then printf "\t\t\t\t\t\t\t\n" >> ${outdir}/${region}.${tppfx}_${testpfx}.${varianttype}.stats
else grep -v "Sample" ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.stats >> ${outdir}/${region}.${tppfx}_${testpfx}.${varianttype}.stats
fi
done < ${samplelist}
done < ${regionconfig}
paste ${outdir}/*.${tppfx}_${testpfx}.${varianttype}.stats > ${outdir}/total_stats.${tppfx}_${testpfx}.${varianttype}


# boxplot data
printf "Sample\tValue\tType\tScope\tMethod\n" >> ${outdir}/total_stats.${tppfx}_${testpfx}.${varianttype}_boxplot
method=${tppfx}"_vs_"${testpfx}
while read sample tpvcf testvcf
do
while read region bed
do
printf ${sample}"\t"`grep -v "Sample" ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.stats|cut -f5`"\tPrecision\t"${region}"\t"${method}"\n" >> ${outdir}/total_stats.${tppfx}_${testpfx}.${varianttype}_boxplot
printf ${sample}"\t"`grep -v "Sample" ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.stats|cut -f6`"\tRecall\t"${region}"\t"${method}"\n" >> ${outdir}/total_stats.${tppfx}_${testpfx}.${varianttype}_boxplot
printf ${sample}"\t"`grep -v "Sample" ${outdir}/${region}_${varianttype}/${sample}.${tppfx}_${testpfx}.${region}_${varianttype}.stats|cut -f7`"\tFscore\t"${region}"\t"${method}"\n" >> ${outdir}/total_stats.${tppfx}_${testpfx}.${varianttype}_boxplot
done < ${regionconfig}
done < ${samplelist}





