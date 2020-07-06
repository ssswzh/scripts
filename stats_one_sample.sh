#!/usr/bin/sh
# @Author: Siwen Zhang, zhangsiwen@grandomics.com, May 31 2019 
# History:
#	20190924, change region to bed file


shellname=$0
# Arguments passing
TEMP=`getopt -o b:r:o:h --long bam:,bed:,outdir:,help \
     -n ${shellname} -- "$@"`

if [ $? != 0 ] ; then printf -- '\033[31mError, please try again with correct arguments.\033[0m\n\n'; printf "Usage:\nsh "${shellname}" --bam bam --bed capture.bed --outdir output_dir \n\n" >&2 ; exit 1 ; fi

#set argument order
eval set -- "$TEMP"
while true ; do
	case "$1" in
		-b|--bam) 
			bam="$2"; shift 2;;
		-r|--bed) 
			bed="$2"; shift 2 ;;
		-o|--outdir)
			outdir="$2"; shift 2 ;;
		-h|--help)
			printf -- 'Usage:\nsh '${shellname}' --bam bam --bed capture.bed --outdir output_dir \n\n';
			printf -- '-b|--bam\t\tbam file path\n';
			printf -- '-r|--bed\t\tcapture region bed\n';
			printf -- '-o|--outdir\t\toutput path\n';
			printf -- '-h|--help\n\t\tprint this help message\n';
			printf -- '@Author: Siwen Zhang, zhangsiwen@grandomics.com, May 31 2019 \n\n';
			exit 0;; #shift ; break; exit 1;;
		--) 
			shift ; break; exit 1;;
		?)
			printf -- 'Usage: \n';
			printf -- 'sh '${shellname}' --bam bam --bed capture.bed --outdir output_dir \n\n';
			exit 0;; #shift ; break; exit 1;;
		*) 
			printf -- 'Usage: ';
			printf -- 'sh '${shellname}' --bam bam --bed capture.bed --outdir output_dir \n\n'
			printf -- '@Author: Siwen Zhang, zhangsiwen@grandomics.com, May 31 2019 \n\n';
			exit 0;; #exit 1;;
	esac
done

set -e 

if [ -n ${bam} ] && [ -n ${bed} ] && [ -n ${outdir} ]
then continue
else exit 1
fi

#chrom=`echo ${region}|cut -d':' -f1`
#startpos=`echo ${region}|cut -d':' -f2|cut -d'-' -f1`
#endpos=`echo ${region}|cut -d':' -f2|cut -d'-' -f2`


# Reads and Mapping ratio
sampleID=`basename ${bam}|cut -d'.' -f1`
[ ! -d ${outdir} ] && mkdir ${outdir}
[ -f ${outdir}/${sampleID}.mapping_ratio.txt ] && rm ${outdir}/${sampleID}.mapping_ratio.txt 
[ -f ${outdir}/${sampleID}.coverage.txt ] && rm ${outdir}/${sampleID}.coverage.txt

printf "SampleID\tFastq_pass_reads\tFastq_pass_base\tMapping_wgs_reads\tMapping_ratio\tWgs_bases\tWgs_avg_depth\tTarget_reads_num\tTarget_reads_ratio\tTarget_bases\tTarget_bases_ratio\tTarget_avg_depth\n" >> ${outdir}/${sampleID}.mapping_ratio.txt

barcoded_reads=`samtools view -F 3840 ${bam}|wc -l`
base=`samtools view -F 3840 ${bam}|cut -f10|wc -c `
barcoded_base=$((${base}-${barcoded_reads}))
mean_length=`Rscript -e "${barcoded_base}/${barcoded_reads}"|awk '{print $2}'`
mapping_wgs_reads=`samtools view -F 3844 ${bam}|wc -l`
mapping_ratio=`Rscript -e "${mapping_wgs_reads}/${barcoded_reads}"|awk '{print $2}'`
wgs_bases=`awk 'BEGIN{FS="\t";t=0}{t=t+$3}END{print t}' ${bam}.depth`
wgs_avg_depth=`Rscript -e "${wgs_bases}/300000000"|awk '{print $2}'`
target_reads_num=`samtools view -F 3844 ${bam} -L ${bed}|wc -l`
target_reads_ratio=`Rscript -e "${target_reads_num}/${mapping_wgs_reads}"|awk '{print $2}'`
target_region=`awk 'BEGIN{s=0;e=0}{s=s+$2;e=e+$3+1}END{print e-s}' ${bed}`
target_bases=$((0))
while read chrom start end
do tmp=`awk -v chr=${chrom} -v start=${start} -v end=${end} 'BEGIN{FS="\t";OFS="\t";total=0}{if($1==chr && $2>=start && $2<=end){total=total+$3}}END{print total}' ${bam}.depth`
target_bases=$(($((${target_base}))+$((${tmp}))))
done < ${bed}
#target_bases=`awk -v chr=${chrom} -v start=${startpos} -v end=${endpos} 'BEGIN{FS="\t";OFS="\t";total=0}{if($1==chr && $2>=start && $2<=end){total=total+$3}}END{print total}' ${bam}.depth`
target_bases_ratio=`Rscript -e "${target_bases}/${wgs_bases}"|awk '{print $2}'`
target_avg_depth=`Rscript -e "${target_bases}/${target_region}"|awk '{print $2}'`

printf $sampleID"\t"$barcoded_reads"\t"$barcoded_base"\t"$mapping_wgs_reads"\t"$mapping_ratio"\t"$wgs_bases"\t"$wgs_avg_depth"\t"$target_reads_num"\t"$target_reads_ratio"\t"$target_bases"\t"$target_bases_ratio"\t"$target_avg_depth"\n" >> ${outdir}/${sampleID}.mapping_ratio.txt


# Coverage
printf "Sample\tAvg_depth\tMedian\tCoverage%%\tCov_4x%%\tCov_10x%%\tCov_30x%%\tCov_100x%%\n" >> ${outdir}/${sampleID}.coverage.txt

sampleID=`basename ${bam}|cut -d'.' -f1`
if [ ! -d ${outdir}/${sampleID} ]
then mkdir ${outdir}/${sampleID}
fi
printf ${chrom}"\t"${startpos}"\t"${endpos} > ${outdir}/${sampleID}/capture.bed
/export/software/conda/miniconda3/bin/bamdst -p ${outdir}/${sampleID}/capture.bed -o ${outdir}/${sampleID} ${bam}
gunzip -c ${outdir}/${sampleID}/region.tsv.gz > ${outdir}/${sampleID}/region.tsv
printf ${sampleID}"\t"`tail -1 ${outdir}/${sampleID}/region.tsv |cut -f4`"\t"`tail -1 ${outdir}/${sampleID}/region.tsv |cut -f5`"\t" >> ${outdir}/${sampleID}.coverage.txt
less ${outdir}/${sampleID}/coverage.report |grep "\[Target\]\ Coverage"|awk '{print $4}'|tr '\n' '\t' >> ${outdir}/${sampleID}.coverage.txt
printf "\n" >> ${outdir}/${sampleID}.coverage.txt

