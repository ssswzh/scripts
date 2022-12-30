#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen
# History:
#   20200221, manuscript

import argparse
from Bio import SeqIO


def GetArgs():
    parser = argparse.ArgumentParser(description='Script for replace FASTA with vcf, trim fasta and print out sequence in fixed length')
    parser.add_argument('--vcf', dest='vcf', help='vcf file', required=True) 
    parser.add_argument('--ref', dest='fa', help='FASTA file', required=True)
    parser.add_argument('--out', dest='out', help='output file name', required=True)
    parser.add_argument('--start', dest='start', help='trim bases before start', required=True) 
    parser.add_argument('--end', dest='end', help='trim bases after end', required=True) 
    parser.add_argument('--prefix', dest='prefix', help='sequence ID written in fasta', required=True) 
    args = parser.parse_args()
    return args


def openVCF(vcffile):
    result = {}
    vcf = open(vcffile)
    for line in vcf:
        ele = line.strip().split("\t")
        if "#" in ele[0]:
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = ele
        if CHROM not in result.keys():
            result[CHROM] = {}
        result[CHROM][POS] = ALT
    vcf.close()
    return result


def chunkstring(string, length):
    result = (string[0+i:length+i] for i in range(0, len(string), length))
    return result


def processFASTAwithVCF(vcffile, fastafile, outfile, prefix, start, end):
    vcf = openVCF(vcffile)
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    out = open(outfile, "w")
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seq_list = list(sequence)
        for pos in vcf[name].keys():
            seq_list[int(pos)-1] = vcf[name][pos]
        new_sequence = ''.join(seq_list)
        new_sequence = new_sequence[int(start):int(end)]
        result = '\n'.join(list(chunkstring(new_sequence, 70)))
       	out.write(">"+prefix+"\n")
        out.write(result+"\n")
    out.close()

def main():
    args = GetArgs()
    vcf, ref, out, prefix, start, end = args.vcf, args.fa, args.out, args.prefix, args.start, args.end
    processFASTAwithVCF(vcf, ref, out, prefix, start, end)


if __name__ == '__main__':
    main()
