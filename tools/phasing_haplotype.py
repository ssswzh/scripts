#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen
# History:
#     20191212, first release

import argparse

def GetArgs():
    parser = argparse.ArgumentParser(description='Script for obtaining phased SNPs in vcf and output two haplotype SNP result.')
    parser.add_argument('--vcf', dest='vcf', help='phased vcf file', required=True) 
    parser.add_argument('--homo', dest='homo', default=True, help='print homozygosity SNPs [True|False], default=True') 
    parser.add_argument('--out', dest='out', help='output file name', required=True)
    args = parser.parse_args()
    return args

def Phase2Haplotype(vcf, prefix, homo):
    vcffile = open(vcf)
    outfile = open(prefix,"w")
    for line in vcffile:
        if not line.startswith("#"):
            ele = line.strip().split("\t")
            pos = str(ele[0]+":"+ele[1])
            ref = ele[3]
            alt = ele[4]
            geno = ele[9].split(":")[0]
            if "|" in geno:
                outfile.write(pos+"\t")
                if geno.split("|")[0] == "0":
                    outfile.write(ref+"\t")
                    outfile.write(alt+"\n")
                else:
                    outfile.write(alt+"\t")
                    outfile.write(ref+"\n")
            elif "/" in geno:
                if homo:
                    geno0 = geno.split("/")[0]
                    geno1 = geno.split("/")[1]
                    if geno0=="1" and geno1=="1":
                        outfile.write(pos+"\t"+alt+"\t"+alt+"\n")
                    
    vcffile.close()
    outfile.close()


def main():
    args = GetArgs()
    Phase2Haplotype(args.vcf, args.out, args.homo)

if __name__=="__main__":
    main()
