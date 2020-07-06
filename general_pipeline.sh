#!/usr/bin/bash

# Created by Zhang Siwen
# History:
#       20190425, first version pipeline
#	20190815, change region to bed, change clairvoyante command
#	20190924, change clair command back to single region SNP calling, remove parallel running
#       20191018, run clair command parallel and merge result
#       20191209, add thread, create folder for clair result
#       20191216, add medaka_variant command for SNP calling
#       20200108, allow custom reference fasta file


if [ $# -lt 7 -o $# -gt 8 ];
    then
    printf "Usage: sh $0 <raw fastq> <outdir> <output prefix> <capture bed> <platform> <ref> <thread> <optional medaka_model> > run.sh\n\n";
    printf "Arguments: \n"
    printf "raw fastq: fastq or fasta file \n"
    printf "outdir: output directory, both absolute path and relative path \n"
    printf "output prefix: output prefix of result file \n"
    printf "capture bed: capture region bed file \n"
    printf "platform: [pb], [ont] \n"
    printf "ref: [hg19], [hg38]. Or any fasta format file \n"
    printf "thread: number of thread \n"
    printf "medaka_model: [r941_min_fast], [r941_min_high], [r941_prom_fast], [r941_prom_high], [r10_min_high], [r941_min_diploid_snp] \n\n"
    exit 1;
fi

# Arguments

sample=$1
#ref_version=$2
outdir=$2
outprefix=$3
bed=$4 
platform=$5
ref_version=$6
thread=$7
medaka_model=$8

annovarDB="/export/database/annovar/humandb"

if [ ${ref_version} == "hg19" ]
	then ref="/export/database/pub_database/hg19/human_g1k_v37.fasta"; chrom_len="/export/home/zhangsiwen/database/hg19/GRCh37.chrom.txt";
        printf "printf \"Choose hg19 as reference, /export/database/pub_database/hg19/human_g1k_v37.fasta \\\n\" \n"
elif [ ${ref_version} == "hg38" ]
	then ref="/export/home/zhangsiwen/database/hg38/hg38.fa";chrom_len="/export/home/zhangsiwen/database/hg38/GRCh38.chrlen.txt";
        printf "printf \"Choose hg38 as reference, /export/home/zhangsiwen/database/hg38/hg38.fa \\\n\" \n"
elif [[ ${ref_version} =~ \.fa$ ]] || [[ ${ref_version} =~ \.fasta$ ]]
        then if [ -f ${ref_version} ]
                then ref=${ref_version};chrom_len="negative"; echo "printf \"Choose '${ref_version}' as reference. \\\n\""
                else echo "Error: File '${ref_version}' not found."; exit 1;
             fi 
else echo "Error: reference version error."; echo "Please choose from [hg19] and [hg38], or provide a sequence file ends with .fa or .fasta."; exit 1;
fi


printf "

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=\$current_command; current_command=\$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo \"\\\\\"\${last_command}\\\" command filed with exit code \$?.\\\"' EXIT
\n"


# Alignment

printf "if [ ! -d ${outdir} ]
	then mkdir -p ${outdir};
fi \n\n"

printf "/export/software/conda/miniconda3/bin/minimap2 --MD -L -Y -t ${thread} --secondary=no -ax map-${platform} ${ref} ${sample} |samtools view -b - |samtools sort -m 4G -o ${outdir}/${outprefix}.bam - && samtools index ${outdir}/${outprefix}.bam \n"
printf "samtools depth -d 0 ${outdir}/${outprefix}.bam > ${outdir}/${outprefix}.bam.depth \n"
printf "samtools view -h -b ${outdir}/${outprefix}.bam -L ${bed} > ${outdir}/${outprefix}.capture.bam \n"
printf "samtools index ${outdir}/${outprefix}.capture.bam \n"
printf "\n\n"


# Statistics
#printf "if [ -d ${outdir}/stats ]; then rm -rf ${outdir}/stats; fi \n"
#printf "mkdir ${outdir}/stats \n"
#printf "sh /export/home/zhangsiwen/scripts/pipeline/stats_one_sample.sh --bam ${outdir}/${outprefix}.bam --bed ${bed} --outdir ${outdir}/stats \n"
printf "sh /export/home/zhangsiwen/scripts/pipeline/Capture_Evaluation_json.sh ${outdir}/${outprefix}.bam ${outdir}/Capture_Evaluation ${outprefix}_capture ${platform} ${ref_version} ${thread} ${bed} > ${outdir}/capture_evaluation.json \n"
printf "printf \"java -Dconfig.file=/export/home/jinhongshuai/capture_WDL/BJ_huaweicloud_SGE.conf -jar /export/software/cromwell-42/cromwell-42.jar run --inputs ${outdir}/capture_evaluation.json /export/home/jinhongshuai/capture_WDL/Capture_Evaluation_v0.1.wdl \" >> ${outdir}/run_capture_evaluation.sh \n"
printf "sh ${outdir}/run_capture_evaluation.sh \n"
printf "\n\n"


# plots
printf "unset PYTHONPATH
python3 /export/home/zhangsiwen/scripts/plot/draw_capture_depth.py --indepth ${outdir}/${outprefix}.bam.depth --bed ${bed} --outfile ${outdir}/${outprefix}.capture_depth \n"
if [ ${chrom_len} != "negative" ]
then 
printf "Rscript /export/home/zhangsiwen/scripts/plot/draw_wgs_depth_distribution.r ${outdir}/${outprefix}.bam.depth ${chrom_len} ${outdir}/${outprefix}.wgs.pdf \n"
fi
printf "\n\n"



# SV calling

printf "/export/home/hanyue/software/Sniffles/bin/sniffles-core-1.0.11/sniffles -t ${thread} -m ${outdir}/${outprefix}.bam -v ${outdir}/${outprefix}.sniffles.vcf --report_BND -n -1 #--ignore_sd -s 1 --genotype -q 0 \n\n"


# Annotation
printf "unset PYTHONPATH \n"
printf "
python3 /export/home/zhangsiwen/scripts/SV_Anno/sv_filter.py --vcf ${outdir}/${outprefix}.sniffles.vcf --min_re 10 --prefix ${outprefix} --outdir ${outdir}
python3 /export/home/zhangsiwen/scripts/SV_Anno/svannot.py --bed ${outdir}/${outprefix}.re10.bed --prefix ${outprefix}.re10 --db ${annovarDB} --annovar /export/software/annovar/table_annovar.pl --outdir ${outdir} 
python3 /export/home/zhangsiwen/scripts/SV_Anno/svan.py --vcf ${outdir}/${outprefix}.re10.bed --prefix ${outprefix}.re10 --outdir ${outdir} --pub_dbdir /export/database/public_sv_database/v0.0.1 --local_dbdir /export/database/grand_sv_database/ONT/database_v0.0.1 
python3 /export/home/zhangsiwen/scripts/SV_Anno/table_maker.py ${outdir}/${outprefix}.re10.hg19.omim_chpo.xls ${outdir}/${outprefix}.re10.svan.tsv ${outdir}/${outprefix}.re10.hg19.final.svannot.xls
python3 /export/home/zhangsiwen/scripts/SV_Anno/tab2excel.py --input_file ${outdir}/${outprefix}.re10.hg19.final.svannot.xls --output_file ${outdir}/${outprefix}.re10.hg19.final.svannot.xlsx \n"
printf "\n\n"


# SNP calling

# 1. bcftools
printf "#bcftools mpileup -Ou -T ${bed} -f ${ref} ${outdir}/${outprefix}.bam |bcftools call -Ov -mv -T ${bed} > ${outdir}/${outprefix}.bcftools.vcf & \n\n"

# 2. freebayes
printf "#freebayes -f ${ref} ${outdir}/${outprefix}.bam -t ${bed} > ${outdir}/${outprefix}.freebayes.vcf \n\n"


# 4. medaka
if [ -z ${medaka_model} ]
then printf "# Not using medaka_variant caller"
else printf "while read chrom start end
do
printf \"medaka_variant -i ${outdir}/${outprefix}.bam -f ${ref} -r \${chrom}:\${start}-\${end} -o ${outdir}/${outprefix}_medaka_\${chrom}_\${start}_\${end} -s ${medaka_model} -m ${medaka_model} -T 0.04 -t ${thread} -b 150 \\\n\" >> ${outdir}/run_${outprefix}.medaka.sh
done < ${bed}
cat ${outdir}/run_${outprefix}.medaka.sh|/home/jinhongshuai/.conda/envs/parallel.env/bin/parallel -j 4
vcfcat ${outdir}/${outprefix}_medaka_*/round_1_phased.vcf| vcfstreamsort > ${outdir}/${outprefix}.medaka_merge.vcf
\n"
printf "/usr/bin/perl /export/software/annovar/table_annovar.pl ${outdir}/${outprefix}.medaka_merge.vcf ${annovarDB} -buildver hg19 -out ${outdir}/${outprefix}.medaka_merge.annovar -otherinfo -nastring . -remove -protocol refGene,wgEncodeGencodeBasicV19,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,rmsk,wgEncodeBroadHmmGm12878HMM,wgEncodeBroadHmmH1hescHMM,wgEncodeBroadHmmHepg2HMM,wgEncodeBroadHmmHuvecHMM,wgEncodeBroadHmmK562HMM,gwasCatalog,phastConsElements46way,avsnp147,clinvar_20190211,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,gnomad_exome,esp6500siv2_all,exac03,dbnsfp33a,gerp++gt2,cg69,intervar_20170202 --operation g,g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --argument '-hgvs  --splicing_threshold 20 ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' -vcfinput \n\n"

fi


# 3. clairvoyante
if [ ${platform} == "ont" ]
then printf "model=fullv3-ont-ngmlr-hg001-hg19/learningRate1e-3.epoch999; threshold=0.25 \n"
else printf "model=fullv3-pacbio-ngmlr-hg001+hg002+hg003+hg004-hg19/learningRate1e-3.epoch100; threshold=0.2 \n"
fi
printf "mkdir ${outdir}/${outprefix}_clair\n"
printf "source /export/software/conda/miniconda3/envs/Clairvoyante/Clairvoyante.bashrc
while read chrom start end
do
printf \"/usr/bin/python /export/software/conda/miniconda3/envs/Clairvoyante/bin/clairvoyante/callVarBam.py --chkpnt_fn /export/software/conda/miniconda3/envs/Clairvoyante/trainedModels/\${model} --bam_fn ${outdir}/${outprefix}.bam --ref_fn ${ref} --sampleName ${outprefix} --call_fn ${outdir}/${outprefix}_clair/${outprefix}.clair.\${chrom}_\${start}_\${end}.vcf --ctgName \${chrom} --ctgStart \${start} --ctgEnd \${end} --threads ${thread} --threshold \${threshold} --minCoverage 4 \\\n\" >> ${outdir}/run_${outprefix}.clair.sh
done < ${bed}

cat ${outdir}/run_${outprefix}.clair.sh|/home/jinhongshuai/.conda/envs/parallel.env/bin/parallel -j 8
vcfcat ${outdir}/${outprefix}_clair/${outprefix}.clair.*.vcf | vcfstreamsort > ${outdir}/${outprefix}.clair_merge.vcf #| bgziptabix ${outdir}/${outprefix}.clair_merge.vcf.gz
\n\n"

#printf "/usr/bin/perl /export/software/annovar/table_annovar.pl ${outdir}/${outprefix}.bcftools.vcf /export/database/annovar/humandb/ -buildver hg19 -out ${outdir}/${outprefix}.bcftools.annovar -otherinfo -nastring . -remove -protocol refGene,wgEncodeGencodeBasicV19,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,rmsk,wgEncodeBroadHmmGm12878HMM,wgEncodeBroadHmmH1hescHMM,wgEncodeBroadHmmHepg2HMM,wgEncodeBroadHmmHuvecHMM,wgEncodeBroadHmmK562HMM,gwasCatalog,phastConsElements46way,avsnp147,clinvar_20190211,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,gnomad_exome,esp6500siv2_all,exac03,dbnsfp33a,gerp++gt2,cg69,intervar_20170202 --operation g,g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --argument '-hgvs  --splicing_threshold 20 ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' -vcfinput \n\n"

printf "/usr/bin/perl /export/software/annovar/table_annovar.pl ${outdir}/${outprefix}.clair_merge.vcf ${annovarDB} -buildver hg19 -out ${outdir}/${outprefix}.clair_merge.annovar -otherinfo -nastring . -remove -protocol refGene,wgEncodeGencodeBasicV19,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,rmsk,wgEncodeBroadHmmGm12878HMM,wgEncodeBroadHmmH1hescHMM,wgEncodeBroadHmmHepg2HMM,wgEncodeBroadHmmHuvecHMM,wgEncodeBroadHmmK562HMM,gwasCatalog,phastConsElements46way,avsnp147,clinvar_20190211,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,gnomad_exome,esp6500siv2_all,exac03,dbnsfp33a,gerp++gt2,cg69,intervar_20170202 --operation g,g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --argument '-hgvs  --splicing_threshold 20 ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' -vcfinput 
\n\n"



# Combine results
#printf "#unset PYTHONPATH \n"
#printf "#python /export/home/zhangsiwen/scripts/merge_bcftools_freebayes_cnn_annotation.py --bcftools ${outdir}/${outprefix}.bcftools.vcf --freebayes ${outdir}/${outprefix}.freebayes.vcf --cnn ${outdir}/${outprefix}.clairvoyante.vcf --output ${outdir}/${outprefix}.merged \n\n"



