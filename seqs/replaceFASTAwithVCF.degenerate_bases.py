#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen
# History:
#   20200221, manuscript
#   20221230, add --Nbed
#   20230103, change --vcf and --Nbed to optional, and at least provide one of them
#   20230105, for multi mutations in the same site, calculate degenerated bases
#   20230106, calculate similarity to reference, because blastn will ignore mismatches in two ends of sequence


import argparse
import sys
from Bio import SeqIO


DegeneratedBases = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AC': 'M',
    'AG': 'R',
    'AT': 'W',
    'CG': 'S',
    'CT': 'Y',
    'GT': 'K',
    'ACG': 'V',
    'ACT': 'H',
    'AGT': 'D',
    'CGT': 'B',
    'ACGT': 'N',
}


def GetArgs():
    parser = argparse.ArgumentParser(description='Script for replacing FASTA with vcf OR/AND bed, optional print out sequence in fixed length')
    parser._action_groups.pop()
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--ref', dest='fa', help='FASTA file', required=True)
    optional = parser.add_argument_group('Optional arguments') 
    optional.add_argument('--vcf', dest='vcf', default=None, help='vcf file', required=False) 
    optional.add_argument('--Nbed', dest='Nbed', default=None, help='give a bed to change sequence in the reference to Ns', required=False) 
    optional.add_argument('--out', dest='out', default='new_sequence', help='output file prefix, default %(default)s', required=False)
    optional.add_argument('--width', dest='width', default=None, type=int, help='print sequence every X bases in a row, default %(default)s', required=False) 
    optional.add_argument('--prefix', dest='prefix', default='new_sequence', help='sequence ID written in fasta, default %(default)s', required=False) 
    optional.add_argument('--deg', dest='deg', default=0.2, type=float, help='AF cutoff for degenerated base calculation, only used in VCF file, default %(default)s', required=False) 
    args = parser.parse_args()
    return args


def readVCF(vcffile, degenerate=0.2): # 20230105
    result = {}
    vcf = open(vcffile)
    ref_pos = [] # store sites which should be REF
    for line in vcf:
        ele = line.strip().split("\t")
        if "#" in ele[0]:
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,SAMPLE = ele
        POS = int(POS)
        # already check reference af and stored position
        if POS in ref_pos:
            continue
        if CHROM not in result.keys():
            result[CHROM] = {}
        if POS not in result[CHROM].keys():
            result[CHROM][POS] = []
        detail = dict(zip(FORMAT.strip().split(':'), SAMPLE.strip().split(':')))
        # check reference af
        ref_ad = int(detail['AD'].split(',')[0])
        dp = int(detail['DP'])
        ref_af = float(ref_ad/dp)
        if ref_af >= 0.5:
            del result[CHROM][POS]
            ref_pos.append(POS)
            continue
        elif ref_af >= degenerate:
            if REF not in result[CHROM][POS]:
                result[CHROM][POS].append(REF)
        # check alt af
        alt_af = float(int(detail['VD'])/dp)
        if alt_af >= degenerate:
            result[CHROM][POS].append(ALT)
    # get all sites and mutations
    for chrom in result.keys():
        for pos in result[chrom].keys():
            mut = result[chrom][pos]
            if type(mut) == list:
                mut.sort()
                deg_key = ''.join(mut)
                result[chrom][pos] = DegeneratedBases[deg_key]
    vcf.close()
    return result


def readBed(Nbed): # 20221230
    bedfile = open(Nbed)
    sites = {}
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


def MergeSites(vcf_sites, nbed_sites): # 20230106
    # returns merged_site_dict, mutation_site_number, Nbed_site_number
    nbed_site_num = sum([len(nbed_sites[key]) for key in nbed_sites.keys()])
    mut_sites = sum([len(vcf_sites[key]) for key in vcf_sites.keys()])
    if vcf_sites == {}:
        return nbed_sites, mut_sites, nbed_site_num
    if nbed_sites == {}:
        return vcf_sites, mut_sites, nbed_site_num
    merged_sites = nbed_sites.copy()
    mut_sites = 0
    for chrom in vcf_sites.keys():
        if chrom in merged_sites.keys():
            for pos in vcf_sites[chrom].keys():
                if pos not in merged_sites[chrom].keys():
                    merged_sites[chrom][pos] = vcf_sites[chrom][pos]
                    mut_sites += 1
        else:
            merged_sites[chrom] = vcf_sites[chrom]
            mut_sites = mut_sites + len(vcf_sites[chrom])
    return merged_sites, mut_sites, nbed_site_num


def chunkstring(string, length):
    result = (string[0+i:length+i] for i in range(0, len(string), length))
    return result


def processFASTAwithVCF(sites, mut_sites, nbed_site_num, fastafile, outfile, prefix, width):
    fasta_sequences = SeqIO.parse(open(fastafile),'fasta')
    out = open(outfile+'.fa', "w")
    total_len = 0 # 20230106
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        total_len += len(sequence) # 20230106
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
    # similarity calculation 20230106
    similarity = 100-float(100*mut_sites/(total_len-nbed_site_num))
    out_similarity = open(outfile+'.similarity', "w")
    out_similarity.write('SequenceLength\tGapPosition\tMutationPosition\tSimilarity\n')
    out_similarity.write('%s\t%s\t%s\t%.4f%%\n' % (total_len, nbed_site_num, mut_sites, similarity))
    out_similarity.close()
    

def main():
    args = GetArgs()
    if args.vcf == None and args.Nbed == None:
        sys.exit("Please provide either --vcf or --Nbed, or both of them.")
    if args.vcf != None: # 20230103
        vcf_sites = readVCF(args.vcf, args.deg)
        print(vcf_sites)
    else:
        vcf_sites = {}
    if args.Nbed != None:
        nbed_sites = readBed(args.Nbed)
        print(nbed_sites)
    else:
        nbed_sites = {}
    merged_sites, mut_site_num, nbed_site_num = MergeSites(vcf_sites, nbed_sites) # 20230106
    out = args.out.strip().replace('.fa','').replace('.fasta','') if args.out.strip().endswith('.fa') or args.out.strip().endswith('.fasta') else args.out
    processFASTAwithVCF(merged_sites, mut_site_num, nbed_site_num, args.fa, out, args.prefix, args.width)


if __name__ == '__main__':
    main()
