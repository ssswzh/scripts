#!/bin/bash

# Created by Zhang Siwen, zhangsiwen@grandomics.com
# History:
#     20200811, manuscript
#     20200812, add colors for stdout message
#     20200813, add checkpoints for each step
#     202008XX, add stats for ccs np, rq
#     20200911, add stats for read length of ccs and individual sample
#     20200914, filter annovar result, limited to exonic and splicing region
#     20201203, add spliceAI, remove several annovar annotation database
# To do:
#     Nothing so far

shellname=$0
script=`which ${shellname}`

# default arguments
platform=pb
reference=hg19
thread=4
outdir=`pwd`
depth_threshold="1,2,3,4,5,6,7,8,9,10,15,20,30,50,100" # to check coverage under these depth

# config
annovarDIR="/sfs-grand-med-clinical/software/annovar/"
svannoDIR="/sfs-grand-med-clinical/software/SV_Anno/"
scriptsDIR="/home/zhangsw/scripts"
annovarDB="/sfs-grand-med-clinical/database/annovar/humandb/" #/export/database/annovar/humandb
pubsvDB="/sfs-grand-med-clinical/database/pubsvDB_v1_20191121/" #/export/database/public_sv_database/v0.0.1
dbSV="/sfs-grand-med-research/database/dbSV/hg19/beds/sniffles/"
clairMODEL="/sfs-grand-med-clinical/database/clair_model/"

# print color
RED='\033[0;31m'
YLW='\033[0;33m'
NC='\033[0m'
bold=`tput bold`


# Arguments passing
TEMP=`getopt -o p:i:b:r:t:T:o:h --long platform:,input:,barcode:,reference:,target:,thread:,outdir:,help \
     -n ${shellname} -- "$@"`

# no arguments detected
if [ $? != 0 ] || [ $# == 0 ] ; 
    then 
        printf "${RED}Error: please try again with correct arguments.${NC}\n\n"; 
        printf "Usage:\nsh "${shellname}" -p platform -i input -b barcode.fa -r ref.fa -t target.bed -T thread -o outdir \n\n" >&2 ; 
        printf "To see detailed arguments, use '-h' or '--help' \n\n";
        exit 1 ; 
fi

# set arguments
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -p|--platform) 
            platform="$2"; shift 2;;
        -i|--input) 
            input="$2"; shift 2;;
        -b|--barcode) 
            barcode="$2"; shift 2;;
        -r|--reference) 
            reference="$2"; shift 2 ;;
        -t|--target) 
            target="$2"; shift 2 ;;
        -T|--thread) 
            thread="$2"; shift 2 ;;
        -o|--outdir)
            outdir="$2"; shift 2 ;;
        -h|--help)
            printf -- 'Usage:\nsh '${shellname}' -p platform -i input -b barcode.fa -r ref.fa -t target.bed -T thread -o outdir \n\n';
            printf -- '-p|--platform\t\t choose from [ont, pb], default [pb] \n';
            printf -- '-i|--input\t\t give a fastq_pass dir or a merged fastq for ont, or subreads.bam for pb \n';
            printf -- '-b|--barcode\t\t barcode fasta file, >P5 should provide for PB \n';
            printf -- '-r|--reference\t\t choose from [hg19, hg38], or give a reference ends with .fasta or .fa, default [hg19] \n';
            printf -- '-t|--target\t\t target bed file \n';
            printf -- '-T|--thread\t\t threads for each sample, default [4] \n\t\t\t\t Four samples will be analysed simultaneously, watch out for memory and threads \n';
            printf -- '-o|--outdir\t\t output path, default current folder \n';
            printf -- '-h|--help\t\t print this help message\n';
            printf -- '@Author: Siwen Zhang, zhangsiwen@grandomics.com, 20200811 \n\n';
            exit 0;; #shift ; break; exit 1;;
        --) 
            shift ; break; exit 1;;
        ?)
            printf -- 'Usage: \n';
            printf -- 'sh '${shellname}' -p platform -i input -b barcode.fa -r ref.fa -t target.bed -T thread -o outdir \n\n';
            printf -- '@Author: Siwen Zhang, zhangsiwen@grandomics.com, 20200811 \n\n';
            exit 0;; #shift ; break; exit 1;;
        *) 
            printf -- 'Usage: ';
            printf -- 'sh '${shellname}' -p platform -i input -b barcode.fa -r ref.fa -t target.bed -T thread -o outdir \n\n'
            printf -- '@Author: Siwen Zhang, zhangsiwen@grandomics.com, 20200811 \n\n';
            exit 0;; #exit 1;;
    esac
done


# exit when any command fails
#set -e
# keep track of the last executed command
#trap 'last_command=\$current_command; current_command=\$BASH_COMMAND' DEBUG
# echo an error message before exiting
#trap 'echo \"\\\\\"\${last_command}\\\" command filed with exit code \$?.\\\"' EXIT


# print arguments
if [ -n ${input} ] && [ -n ${barcode} ] && [ -n ${target} ]
    then 
        echo -e "${RED}${bold}[WARNING]${NC}"
        echo -e "${RED}If you want to re-run pipeline, please remove 'done_XXX' file in output directory first.${NC}\n"
        echo -e "${YLW}${bold}Arguments:${NC}"
        echo "Platform: ${platform}"
        echo "Input: ${input}"
        echo "Barcode: ${barcode}"
        echo "Reference: ${reference}"
        echo "Target: ${target}"
        echo -e "Thread: ${thread} ${bold}[for each sample]${NC}"
        echo -e "OutputDIR: ${outdir}\n"
    else 
        exit 1
fi


# full command
echo -e "${YLW}${bold}Full command:${NC}"
echo -e "sh ${script} -p ${platform} -i ${input} -b ${barcode} -r ${reference} -t ${target} -T ${thread} -o ${outdir}\n"
echo -e "${RED}${bold}[WARNING]${NC}"
echo -e "${RED}Four samples will be analysed simultaneously, watch out total memory and threads.${NC}\n"


# reference check
if [ ${reference} == "hg19" ]
    then 
        ref="/sfs-grand-med-clinical/database/hg19/human_g1k_v37.fasta" #"/export/database/pub_database/hg19/human_g1k_v37.fasta"; 
        #chrom_len="/export/home/zhangsiwen/database/hg19/GRCh37.chrom.txt";
elif [ ${reference} == "hg38" ]
    then 
        ref="/sfs-grand-med-research/database/human/GRCh38.p13_GRCh38_mainChrom/GCF_000001405.39_GRCh38.p13_genomic.main_chrom.fa" #"/export/home/zhangsiwen/database/hg38/hg38.fa";
        #chrom_len="/export/home/zhangsiwen/database/hg38/GRCh38.chrlen.txt";
elif [[ ${reference} =~ \.fa$ ]] || [[ ${reference} =~ \.fasta$ ]]
        then 
            if [ -f ${reference} ]
                then 
                    ref=${reference};#chrom_len="negative";
                else 
                    echo "Error: File '${reference}' not found."; exit 1;
            fi 
else 
    echo -e `date`" ${RED}Error: reference version error.${NC}\n"; 
    echo "Please choose from [hg19] and [hg38], or provide a sequence file ends with .fa or .fasta."; 
    exit 1;
fi
echo -e `date`" ${YLW}Choose ${reference} as reference${NC}, path: ${ref} "


# barcode check
if [ ${platform} == "pb" ]
    then
        p5=`grep ">P5\|>p5" ${barcode}`
        if [ -z ${p5} ]
            then
                echo -e `date`" ${RED}${bold}[Error]${NC}"
                echo -e `date`" ${RED}Did not find P5 in barcode file.${NC}"
                exit 1
        fi
fi


# folder check
if [[ ! "${outdir}" == /* ]]
    then 
        outpath=`pwd`"/"${outdir}
        outdir=${outpath}
        echo -e `date`" ${YLW}Received relative path${NC}, ${outdir}"
    else
        echo -e `date`" ${YLW}Received absolute path${NC}, ${outdir}"
fi

if [ ! -d ${outdir}/demultiplex ]
    then
        mkdir -p ${outdir}/demultiplex ${outdir}/mapping ${outdir}/depth_coverage ${outdir}/sniffles ${outdir}/clair ${outdir}/deepvariant
    else
        echo -e `date`" ${RED}${bold}[WARNING]${NC} Output foler exists, overwrite the results."
fi 
echo "Pipeline started at "`date` >> ${outdir}/start_pipeline



# ccs
if [ ${platform} == "pb" ] && [ ! -f ${outdir}/ccs/done_ccs ]
    then
        # check input format
        if [[ ! ${input} =~ \.bam$ ]]
            then 
                echo -e `date`" ${RED}Invalid input.${NC}"
                echo "If you choose pacbio platform [pb], please provide subreads.bam."
            exit 1
        fi
        # check subreads.bam index file
        if [ ! -f ${input}.pbi ]
            then pbindex ${input}
        fi
        # ccs
        mkdir -p ${outdir}/ccs
        if [ ! -f ${outdir}/ccs/ccs.8.bam ]
            then
                echo -e `date`" ${YLW}Start CCS.${NC}"
                for i in {1..8}
                do
                    echo "ccs --min-passes 5 ${input} ${outdir}/ccs/ccs.${i}.bam --chunk ${i}/8 -j 3 --report-file ${outdir}/ccs/ccs.${i}.report --log-file ${outdir}/ccs/ccs.${i}.log"
                done > ${outdir}/ccs/run.sh
                cat ${outdir}/ccs/run.sh|parallel -j 8
                echo -e `date`" ${YLW}CCS Finished.${NC}"
            else
                echo -e `date`" CCS results exist, ${YLW}skipped.${NC}"
        fi
        # merge ccs bam
        if [ ! -f ${outdir}/ccs/ccs.merged.bam ]
            then
                echo -e `date`" ${YLW}Start merge CCS bams.${NC}"
                samtools merge -@ ${thread} ${outdir}/ccs/ccs.merged.bam ${outdir}/ccs/ccs.*.bam
                samtools index ${outdir}/ccs/ccs.merged.bam
                samtools view ${outdir}/ccs/ccs.merged.bam|awk 'BEGIN{FS="\t";OFS="\t"}{split($13,np,":"); split($14,rq,":"); print $1,np[3],rq[3],length($10)}' > ${outdir}/ccs/ccs.merged.bam.np_rq_rl # modified 20200911
                Rscript ${scriptsDIR}/plot_stats_summary.r ${outdir}/ccs/ccs.merged.bam.np_rq_rl 2 "Number_Pass"
                Rscript ${scriptsDIR}/plot_stats_summary.r ${outdir}/ccs/ccs.merged.bam.np_rq_rl 3 "Read_Quality"
                Rscript ${scriptsDIR}/plot_stats_summary.r ${outdir}/ccs/ccs.merged.bam.np_rq_rl 4 "Read_Length" # added 20200911
                echo -e `date`" ${YLW}Merging CCS bams finished.${NC}"
            else
                echo -e `date`" Merged CCS results exist, ${YLW}skipped.${NC}"
        fi
        echo "CCS finished at "`date` > ${outdir}/ccs/done_ccs
fi


# demultiplex
if [ ${platform} == "ont" ] && [ ! -f ${outdir}/demultiplex/done_demultiplex ]
    then
        # check input is a folder or one fastq file, and demultiplex
        echo -e `date`" ${YLW}Start demultiplex.${NC}"
        if [ -d ${input} ]
            then
                cat ${input}/*fastq|nanoplexer -b ${barcode} -p ${outdir}/demultiplex -
        elif [[ ${input} =~ \.fq$ ]] || [[ ${input} =~ \.fastq$ ]]
            then
                nanoplexer -b ${barcode} -p ${outdir}/demultiplex ${input}
        else
            echo -e `date`" ${RED}Invalid input.${NC}"
            echo "If you choose ONT platform [ont], please provide fastq_pass dir or a merged fastq file ends with .fastq or .fq."
            exit 1
        fi
        echo -e `date`" ${YLW}Demultiplex finished.${NC}"
        echo "Demultiplex finished at "`date` > ${outdir}/demultiplex/done_demultiplex

elif [ ${platform} == "pb" ] && [ ! -f ${outdir}/demultiplex/done_demultiplex ]
    then
        if [ ! -f ${outdir}/demultiplex/sample_bam_list ]
            then
                echo -e `date`" ${YLW}Start demultiplex.${NC}"
                lima --ccs --different --split-bam-named -j ${thread} ${outdir}/ccs/ccs.merged.bam ${barcode} ${outdir}/demultiplex/demultiplex.bam --log-file ${outdir}/demultiplex/demultiplex.log
                if [ `grep "^>" ${barcode}|head -1` == ">P5" -o `grep "^>" ${barcode}|head -1` == ">p5" ]
                    then 
                        ls ${outdir}/demultiplex/demultiplex.*[pP]5*.bam|awk '{split($1,a,"/"); split(a[length(a)],b,"."); split(b[2],c,"-"); print c[3]"\t"a[length(a)]}' > ${outdir}/demultiplex/sample_bam_list
                elif [ `grep "^>" ${barcode}|tail -1` == ">P5" -o `grep "^>" ${barcode}|tail -1` == ">p5" ]
                    then
                        ls ${outdir}/demultiplex/demultiplex.*[pP]5*.bam|awk '{split($1,a,"/"); split(a[length(a)],b,"."); split(b[2],c,"-"); print c[1]"\t"$1}' > ${outdir}/demultiplex/sample_bam_list
                fi
                echo -e `date`" ${YLW}Demultiplex finished.${NC}"
            else
                echo -e `date`" Demultiplex bam results exist, ${YLW}skipped.${NC}"
        fi
        # bam2fastq
        if [ ! -f "${outdir}/demultiplex/`awk 'NR==1{print $1}' ${outdir}/demultiplex/sample_bam_list`.fastq" ]
            then
                echo -e `date`" ${YLW}Start bam2fastq format transform.${NC}"
                cat ${outdir}/demultiplex/sample_bam_list|while read prefix bam
                do
                    bam2fastq ${bam} -u -o ${outdir}/demultiplex/${prefix}
                    awk 'NR%4==2{print length($0)}' ${outdir}/demultiplex/${prefix}.fastq > ${outdir}/demultiplex/${prefix}.fastq.rl # added 20200911
                    Rscript ${scriptsDIR}/plot_stats_summary.r ${outdir}/demultiplex/${prefix}.fastq.rl 1 "Read_Length" # added 20200911
                done
                echo -e `date`" ${YLW}Bam2fastq format transform finished.${NC}"
            else
                echo -e `date`" Demultiplex fastq results exist, ${YLW}skipped.${NC}"
        fi
        echo "Demultiplex finished at "`date` > ${outdir}/demultiplex/done_demultiplex
    else
        echo -e `date`" CCS and demultiplex results exist, ${YLW}skipped.${NC}"
fi



# sample id list
ls ${outdir}/demultiplex/*fastq |grep -v "unclassified\|none"|awk '{split($1,a,"/"); split(a[length(a)],b,"."); print b[1]}' > ${outdir}/sample_id
sampleID=${outdir}/sample_id


# mapping
if [ ! -f ${outdir}/mapping/done_mapping ]
    then
        echo -e `date`" ${YLW}Start mapping.${NC}"
        for sample in `cat ${sampleID}`
        do
            # mapping
            printf "minimap2 --MD -L -Y -t ${thread} --secondary=no -ax map-${platform} ${ref} ${outdir}/demultiplex/${sample}.fastq |samtools view -b - |samtools sort -m 4G -o ${outdir}/mapping/${sample}.bam -  &&  " 
            printf "samtools index ${outdir}/mapping/${sample}.bam  &&  "
            printf "samtools depth -d 0 ${outdir}/mapping/${sample}.bam -b ${target} > ${outdir}/mapping/${sample}.target.depth \n" 
        done > ${outdir}/mapping/run_mapping.sh
        cat ${outdir}/mapping/run_mapping.sh |parallel -j 4
        echo -e `date`" ${YLW}Mapping finished.${NC}"
        echo "Alignment finished at "`date` > ${outdir}/mapping/done_mapping
    else
        echo -e `date`" Alignment results exist, ${YLW}skipped.${NC}"
fi


# stats
if [ ! -f ${outdir}/depth_coverage/done_stats ]
    then
        echo -e `date`" ${YLW}Start calculate mapping stats.${NC}"
        for sample in `cat ${sampleID}`
        do
            printf "mosdepth --by ${target} -t ${thread} -T ${depth_threshold} ${outdir}/depth_coverage/${sample} ${outdir}/mapping/${sample}.bam  &&  samtools stats -@ ${thread} ${outdir}/mapping/${sample}.bam > ${outdir}/mapping/${sample}.bam.stat  &&  "
            printf "python ${scriptsDIR}/FinalStat.py --prefix ${sample} --MapStat ${outdir}/mapping/${sample}.bam.stat --bam ${outdir}/mapping/${sample}.bam --capture_bed ${target} --mosdepth_regions ${outdir}/depth_coverage/${sample}.regions.bed.gz --mosdepth_thresholds ${outdir}/depth_coverage/${sample}.thresholds.bed.gz --mosdepth_per_base ${outdir}/depth_coverage/${sample}.per-base.bed.gz --outdir ${outdir}/depth_coverage  &&  " 
            printf "python3 ${scriptsDIR}/draw_capture_depth.py --indepth ${outdir}/mapping/${sample}.target.depth --bed ${target} --outfile ${outdir}/depth_coverage/${sample}.target.depth \n" 
        done > ${outdir}/depth_coverage/run_stats.sh
        cat ${outdir}/depth_coverage/run_stats.sh |parallel -j 4
        echo -e `date`" ${YLW}Mapping stats calculation finished.${NC}"
        echo "Stats finished at "`date` > ${outdir}/depth_coverage/done_stats
    else
        echo -e `date`" Alignment stats results exist, ${YLW}skipped.${NC}"
fi


# sniffles
if [ ! -f ${outdir}/sniffles/done_sniffles ]
    then
        for sample in `cat ${sampleID}`
        do
            if [ ${platform} == "ont" ]
                then re="10"
                    printf "sniffles -t ${thread} -m ${outdir}/mapping/${sample}.bam -v ${outdir}/sniffles/${sample}.vcf --report_BND -n -1 --min_support ${re} && " 
            elif [ ${platform} == "pb" ]
                then re="5"
                    printf "sniffles -t ${thread} -m ${outdir}/mapping/${sample}.bam -v ${outdir}/sniffles/${sample}.vcf --report_BND -n -1 --min_support ${re} --ccs_reads && "
            fi
            printf "python3 ${svannoDIR}/sv_filter.py --vcf ${outdir}/sniffles/${sample}.vcf --min_re ${re} --prefix ${sample} --outdir ${outdir}/sniffles/ && " 
            printf "python3 ${svannoDIR}/svannot.py --bed ${outdir}/sniffles/${sample}.re${re}.bed --prefix ${sample}.re${re} --db ${annovarDB} --annovar ${annovarDIR}/table_annovar.pl --outdir ${outdir}/sniffles \n" 
            #printf "python3 ${svannoDIR}/svan.py --vcf ${outdir}/sniffles/${sample}.re${re}.bed --prefix ${sample}.re${re} --outdir ${outdir}/sniffles --pub_dbdir ${pubsvDB} --local_dbdir ${dbSV}  && " 
            #printf "python3 ${svannoDIR}/table_maker.py ${outdir}/sniffles/${sample}.re${re}.hg19.omim_chpo.xls ${outdir}/sniffles/${sample}.re${re}.svan.tsv ${outdir}/sniffles/${sample}.re${re}.hg19.final.svannot.xls \n" 
        done > ${outdir}/sniffles/run_sniffles.sh
        cat ${outdir}/sniffles/run_sniffles.sh |parallel -j 4
        echo -e `date`" ${YLW}Sniffles finished.${NC}"
        echo "Sniffles finished at "`date` > ${outdir}/sniffles/done_sniffles
    else
        echo -e `date`" Sniffles results exist, ${YLW}skipped.${NC}"
fi


# SNV calling
# clair
if [ ! -f ${outdir}/clair/done_clair ]
    then
        for sample in `cat ${sampleID}`
        do
            #printf "clair.py callVarBam --chkpnt_fn ${clairMODEL}/${platform}_model --bam_fn ${outdir}/mapping/${sample}.bam --ref_fn ${ref} --call_fn ${outdir}/clair/${sample}.vcf --sampleName ${sample} --pysam_for_all_indel_bases --threshold 0.2 --threads ${thread} --minCoverage 4 --qual 100 --bed_fn ${target} \n"
            cut -f1-3 ${target}|while read chrom start end
            do
                printf "docker run --ipc=host --user=`id -u`:`id -g` -v ${outdir}:/data/ -v `dirname ${ref}`:/database/ -v ${clairMODEL}:/model/ swr.cn-north-1.myhuaweicloud.com/grand-med-clinical/clair:2.1.0 clair.py callVarBam --chkpnt_fn /model/${platform}_model --bam_fn /data/mapping/${sample}.bam --ref_fn /database/`basename ${ref}` --call_fn /data/clair/${sample}.vcf --sampleName ${sample} --pysam_for_all_indel_bases --threshold 0.2 --threads ${thread} --minCoverage 4 --qual 100 --ctgName ${chrom} --ctgStart ${start} --ctgEnd ${end} \n"
            done
        done > ${outdir}/clair/run_clair.sh
        cat ${outdir}/clair/run_clair.sh |parallel -j 4
        echo -e `date`" ${YLW}Clair finished.${NC}"
        echo "Clair finished at "`date` > ${outdir}/clair/done_clair
    else
        echo -e `date`" Clair results exist, ${YLW}skipped.${NC}"
fi


# deepvariant
if [ ! -f ${outdir}/deepvariant/done_deepvariant ]
    then
        if [ ${platform} == "ont" ]
            then
                for sample in `cat ${sampleID}`
                do
                    printf "docker run --ipc=host --user=`id -u`:`id -g` -v ${outdir}:/data/ -v `dirname ${ref}`:/database/ -v `dirname ${target}`:/target/ swr.cn-north-1.myhuaweicloud.com/grand-med-clinical/pepper_deepvariant_cpu:latest /opt/run_pepper_deepvariant.sh -f /database/`basename ${ref}` -b /data/mapping/${sample}.bam -o /data/deepvariant/${sample} -t ${thread} -r /target/`basename ${target}` -x 1  &&  "
                    printf "cp ${outdir}/deepvariant/${sample}/PEPPER_SNP_DEEPVAIRANT_OUPUT.vcf.gz ${outdir}/deepvariant/${sample}.vcf.gz  &&  gunzip ${outdir}/deepvariant/${sample}.vcf.gz \n"
                done > ${outdir}/deepvariant/run_deepvariant.sh
        elif [ ${platform} == "pb" ]
            then
                for sample in `cat ${sampleID}`
                do
                    printf "docker run --ipc=host --user=`id -u`:`id -g` -v ${outdir}:/data/ -v `dirname ${ref}`:/database/ -v `dirname ${target}`:/target/ swr.cn-north-1.myhuaweicloud.com/grand-med-research/deepvariant:0.10.0 /opt/deepvariant/bin/run_deepvariant --model_type=PACBIO --ref=/database/`basename ${ref}` --reads=/data/mapping/${sample}.bam --output_vcf=/data/deepvariant/${sample}.vcf --regions=/target/`basename ${target}` --sample_name=${sample} --num_shards=${thread} \n"
                done > ${outdir}/deepvariant/run_deepvariant.sh
        fi
        cat ${outdir}/deepvariant/run_deepvariant.sh |parallel -j 4
        echo -e `date`" ${YLW}Deepvariant finished.${NC}"
        echo "Deepvariant finished at "`date` > ${outdir}/deepvariant/done_deepvariant
    else
        echo -e `date`" Deepvariant results exist, ${YLW}skipped.${NC}"
fi


# optional medaka
#if [ ${platform} == "ont" ]
#    then
#        if [ ! -f ${outdir}/medaka/done_medaka ]
#            then
#                for sample in `cat ${sampleID}`
#                do
#                    mkdir -p ${outdir}/medaka/${sample}_medaka
#                    cut -f1-3 ${target}|while read chrom start end
#                    do
#                        printf "medaka_variant -i ${outdir}/mapping/${sample}.bam -f ${ref} -r ${chrom}:${start}-${end} -o ${outdir}/medaka/${sample}_medaka/${sample}_medaka_${chrom}_${start}_${end} -s r941_prom_snp_g303 -m r941_prom_variant_g303 -T 0.04 -t ${thread} -b 150 \n" 
#                    done 
#                done > ${outdir}/medaka/run_medaka.sh
#                cat ${outdir}/medaka/run_medaka.sh|parallel -j 4
#                for sample in `cat ${sampleID}`
#                do
#                    vcfcat ${outdir}/medaka/${sample}_medaka_*/round_1.vcf| vcfstreamsort > ${outdir}/medaka/${sample}.vcf
#                done
#                echo -e `date`" ${YLW}Medaka finished.${NC}"
#                echo "Medaka finished at "`date` > ${outdir}/medaka/done_medaka
#            else
#                echo -e `date`" Medaka results exist, ${YLW}skipped.${NC}"
#        fi
#    else
#        echo -e `date`" Medaka not run for pb data."
#fi
#


# spliceAI, added 20201203
if [ ${ref} == "hg19" ]
    then
        spliceref="grch37"
elif [ ${ref} == "hg38" ]
    then
        spliceref="grch38"
fi


# annovar
if [ ! -f ${outdir}/done_annovar ]
    then
        for vcf in `ls ${outdir}/[cd]*/*vcf`
        do
            prefix=`dirname ${vcf}`"/"`basename -s .vcf ${vcf}` # added 20201203
            printf "spliceai -I ${vcf} -O ${prefix}.spliceai.vcf -R ${ref} -A ${spliceref} && " # added 20201203
            #printf "perl ${annovarDIR}/table_annovar.pl ${vcf} ${annovarDB} -buildver ${reference} -out ${vcf} -otherinfo -nastring . -remove -protocol refGene,wgEncodeGencodeBasicV19,cpgIslandExt,cytoBand,wgRna,targetScanS,tfbsConsSites,genomicSuperDups,rmsk,wgEncodeBroadHmmGm12878HMM,wgEncodeBroadHmmH1hescHMM,wgEncodeBroadHmmHepg2HMM,wgEncodeBroadHmmHuvecHMM,wgEncodeBroadHmmK562HMM,gwasCatalog,phastConsElements46way,avsnp147,clinvar_20190211,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,gnomad_exome,esp6500siv2_all,exac03,dbnsfp33a,gerp++gt2,cg69,intervar_20170202 --operation g,g,r,r,r,r,r,r,r,r,r,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f --argument '-hgvs  --splicing_threshold 20 ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,' -vcfinput && awk '\$6==\"exonic\"||\$6==\"splicing\"||\$6==\"Func.refGene\"' ${vcf}.${reference}_multianno.txt > ${vcf}.${reference}_multianno.filter.txt \n" # modified 20201203
            printf "perl ${annovarDIR}/table_annovar.pl ${prefix}.spliceai.vcf ${annovarDB} -buildver ${reference} -out ${prefix}.spliceai -otherinfo -nastring . -remove -protocol refGene,cytoBand,genomicSuperDups,rmsk,clinvar_20190211,1000g2015aug_eas,1000g2015aug_sas,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_all,gnomad_exome,exac03 --operation g,r,r,r,f,f,f,f,f,f,f,f,f --argument '-hgvs  --splicing_threshold 20 ,,,,,,,,,,,,' -vcfinput && awk '\$6==\"exonic\"||\$6==\"splicing\"||\$6==\"Func.refGene\"' ${prefix}.spliceai.${reference}_multianno.txt > ${prefix}.spliceai.${reference}_multianno.filter.txt \n" # modified 20201203
        done > ${outdir}/run_annovar.sh
        cat ${outdir}/run_annovar.sh |parallel -j 4
        echo -e `date`" ${YLW}Annovar finished.${NC}"
        echo "Annotation finished at "`date` > ${outdir}/done_annovar
    else
        echo -e `date`" Annovar results exist, ${YLW}skipped.${NC}"
        echo -e `date`" ${YLW}All steps had generated their results.${NC}"
        echo -e "${RED}${bold}[WARNING]${NC}"
        echo -e "${RED}${bold}Did you forget to remove intermediate results? OR what do you want from me exactly?${NC}"
fi


echo -e `date`" ${YLW}${bold}Pipeline finished.${NC}"
echo "Pipeline finished at "`date` > ${outdir}/done_pipeline

