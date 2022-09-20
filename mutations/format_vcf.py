#!/usr/bin/env python3
# Author: Zhangsiwen
# History:
#     20200708, manuscript, to format variants in vcf files


import argparse
import sys

def GetArgs():
    parser = argparse.ArgumentParser(description='Format VCF files, seperate multiple SNVs in the same record.')
    parser.add_argument('--vcf', dest='vcf', help="vcf file", required=True)
    parser.add_argument('--out', dest='out', help='out vcf file prefix (with or without dir), will generate prefix.format.vcf', required=True) 
    args = parser.parse_args()
    return args


def RegulatePositions(ele, pos, ref, alt, records): 
    if len(ref)==len(alt): 
        if len(ref)==1: # len(ref)=len(alt)=1, ref and alt are both single nucleotide
            records.append('\t'.join(ele))
        else: # len(ref)=len(alt)>1
            for chara in range(len(ref)):
                # base not equal in same position of ref and alt, like AC and GT, keep as single nucleotide variant
                if ref[chara]!=alt[chara]: 
                    ele[1] = str(int(pos)+chara)
                    ele[3] = ref[chara]
                    ele[4] = alt[chara]
                    records.append('\t'.join(ele))
    else:
        ele[4] = alt
        records.append('\t'.join(ele))


def VCFFormatConversion(vcffile, out): 
    vcf = open(vcffile)
    records = []
    for line in vcf:
        if line.startswith("#"): # keep headers
            records.append(line.strip())
        else:
            ele = line.strip().split("\t")
            pos, ref, alt = ele[1], ele[3], ele[4]
            if "," not in alt: # alt variant is single type
                RegulatePositions(ele, pos, ref, alt, records)
            else: # alt variant contains multiple formats, calculate each one of them
                allele = alt.split(",")
                for i in allele:
                    RegulatePositions(ele, pos, ref, i, records)
    outfile = open(out, "w")
    for r in records:
        outfile.write(r+"\n")
    outfile.close()


def main():
    args = GetArgs()
    format_vcf = args.out+".format.vcf"
    VCFFormatConversion(args.vcf, format_vcf)

if __name__=="__main__":
    main()


