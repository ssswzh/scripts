#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen
# History:
#   20200221, manuscript
#   20221230, add --Nbed
#   20230103, change --vcf and --Nbed to optional, and at least provide one of them


import argparse
from Bio import SeqIO


def GetArgs():
    parser = argparse.ArgumentParser(description='Script for replace FASTA with vcf, trim fasta and print out sequence in fixed length')
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--ref', dest='fa', help='FASTA file', required=True)
    optional = parser.add_argument_group('Optional arguments') 
    optional.add_argument('--vcf', dest='vcf', default=None, help='vcf file', required=False) 
    optional.add_argument('--Nbed', dest='Nbed', default=None, help='give a bed to change sequence in the reference to Ns', required=False) 
    optional.add_argument('--out', dest='out', default='new_sequence.fa', help='output file name', required=False)
    optional.add_argument('--width', dest='width', default=None, type=int, help='print sequence every X bases in a row', required=False) 
    optional.add_argument('--prefix', dest='prefix', default='new_sequence', help='sequence ID written in fasta', required=False) 
    args = parser.parse_args()
    return args


def readVCF(vcffile):
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


def readBed(Nbed, sites): # 20221230
    bedfile = open(Nbed)
    for line in bedfile:
        if line.startswith('#'):
            continue
        ele = line.strip().split('\t')
        chrom, start, end = ele[0:3]
        if chrom not in sites.keys():
            sites[chrom] = {}
        for pos in range(int(start)+1, int(end)+1):
            sites[chrom][pos] = 'N'
    return sites


def chunkstring(string, length):
    result = (string[0+i:length+i] for i in range(0, len(string), length))
    return result


def processFASTAwithVCF(sites, fastafile, outfile, prefix, width):
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    out = open(outfile, "w")
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        seq_list = list(sequence)
        for pos in sites[name].keys():
            seq_list[int(pos)-1] = sites[name][pos]
        new_sequence = ''.join(seq_list)
        if width is None:
            result = new_sequence
        else:
            result = '\n'.join(list(chunkstring(new_sequence, width)))
       	out.write(">"+prefix+"\n")
        out.write(result+"\n")
    out.close()

def main():
    args = GetArgs()
    if args.vcf != None: # 20230103
        sites = readVCF(args.vcf)
    else:
        sites = {}
    if args.Nbed != None:
        sites = readBed(args.Nbed, sites)
    processFASTAwithVCF(sites, args.fa, args.out, args.prefix, args.width)


if __name__ == '__main__':
    main()
