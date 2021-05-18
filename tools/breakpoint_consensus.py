#!/usr/bin/env python3
# Author: Zhang Siwen, Han Yue
# History:
#     20210514, first version (lost the original shell script)

import argparse
import pysam
import logging
import os
from decode_mapping_position import ParseBamFile


def GetArgs():
    parser = argparse.ArgumentParser(description='Use SpliceAI model to predict splice site for SVs (only deletions and duplications).')
    parser.add_argument('--prefix', dest='prefix', help='output prefix', required=True) 
    parser.add_argument('--sv', dest='sv', help='sv info, chrom:start-end', required=True) 
    parser.add_argument('--vcf', dest='vcf', help='vcf file', required=True) 
    parser.add_argument('--bam', dest='bam', help='mapping bam file', required=True) 
    parser.add_argument('--out', dest='outpath', default="./", help='output path, current path for default')
    parser.add_argument('--ref', dest='reference', default="hg19", help='reference version, hg19 or hg38, default [hg19]')
    parser.add_argument('--log', dest='log', action='store_true', help='save intermediate results')
    args = parser.parse_args()
    return args


def RefVersion(ref):
    if ref == "hg19":
        return "/sfs-grand-med-research/database/human/humanG1Kv37/human_g1k_v37.fasta"
    elif ref == "hg38":
        return "/sfs-grand-med-research/database/human/GRCh38.p13_GRCh38_mainChrom/GCF_000001405.39_GRCh38.p13_genomic.main_chrom.fa"
    else:
        print("Wrong reference version\n")
        exit()
    

def ExtractReadID(outprefix, chrom, start, end, vcffile, log):
    # outfiles
    if log:
        out_readid = open(outprefix + ".readid", "w")
    # extract read id from vcf
    readid = []
    try:
        vcf = pysam.VariantFile(vcffile)
        vcf.header.filters.add("STRANDBIAS",None,None,"Dummy")
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    for record in vcf:
        if str(chrom) == str(record.chrom) and int(start) == int(record.pos) and int(end) == int(record.stop):
            readid = list(record.info['RNAMES'])
            if log:
                for id in readid:
                    out_readid.write(id + "\n")
                out_readid.close()
    # exit if not find sv
    if readid == []:
        logging.error("Error: " + sv + " not found in " + vcf)
        exit()
    return readid


def ExtractReadBam(outprefix, readid, bamfile):
    # extract reads from bam file
    bam = pysam.Samfile(bamfile, "rb")
    # out fastq file 
    out = open(outprefix + ".fastq", 'w')
    for record in bam:
        if record.qname in readid:
            out.write("@" + record.query_name + "\n")
            out.write(record.seq + "\n")
            out.write("+\n")
            out.write(record.qual + "\n")
            readid.remove(record.qname)
    out.close()


def ConsensusMapping(outprefix, prefix, ref):
    fq = outprefix + ".fastq"
    os.system("spoa {} |sed \"s/Consensus/>{}/\" > {}.consensus.fasta".format(fq, prefix, outprefix))
    os.system("minimap2 --MD -L -Y --secondary=no -ax asm10 {} {}.consensus.fasta > {}.consensus.sam".format(ref, outprefix, outprefix))
    ParseBamFile(outprefix + ".consensus.sam", outprefix + ".consensus.sam.pos", False)


def main():
    args = GetArgs()
    # sv positions
    chrom = args.sv.strip().split(":")[0]
    start = args.sv.strip().split(":")[1].split("-")[0]
    end = args.sv.strip().split(":")[1].split("-")[1]
    outprefix = args.outpath + "/" + args.prefix + "." + chrom + "_" + start + "_" + end
    # get reference
    ref = RefVersion(args.reference)
    # process
    readid = ExtractReadID(outprefix, chrom, start, end, args.vcf, args.log)
    ExtractReadBam(outprefix, readid, args.bam)
    ConsensusMapping(outprefix, args.prefix, ref)
    

if __name__ == '__main__':
    main()

