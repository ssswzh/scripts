#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen
# History:
#   20200225, manuscript

import argparse
from operator import itemgetter
from itertools import groupby


def GetArgs():
    parser = argparse.ArgumentParser(description='Find continuous ranges of position with depth larger than threshold.')
    parser.add_argument('--depth', dest='depth', help='sample.bam.depth', required=True) 
    parser.add_argument('--thre', dest='thre', help='keep pos with threshold >= thre', required=True)
    parser.add_argument('--out', dest='out', help='output file prefix', required=True)
    args = parser.parse_args()
    return args

def ListToRange(alist):
    ranges = []
    for k,g in groupby(enumerate(alist),lambda x:x[0]-x[1]):
        group = (map(itemgetter(1),g))
        group = list(map(int,group))
        ranges.append([group[0],group[-1]])
    return ranges

def FindRanges(infile, cutoff, outfile):
    depth = open(infile)
    pos_keep = {}
    pos_lose = {}
    for line in depth:
        ele = line.strip().split("\t")
        chrom, pos, dp = ele[0], int(ele[1]), int(ele[2])
        if dp >= int(cutoff):
            if chrom not in pos_keep.keys():
                pos_keep[chrom] = []
            pos_keep[chrom].append(pos)
        else :
            if chrom not in pos_lose.keys():
                pos_lose[chrom] = []
            pos_lose[chrom].append(pos)

    out_keep = open(outfile+".above_thre_"+cutoff+".bed","w")
    for key in pos_keep.keys():
        keep_ranges = ListToRange(pos_keep[key])
        for ele in keep_ranges:
            out_keep.write("\t".join([key, str(ele[0]), str(ele[1])]))
            out_keep.write("\n")
    out_keep.close()

    out_lose = open(outfile+".below_thre_"+cutoff+".bed","w")
    for key in pos_lose.keys():
        lose_ranges = ListToRange(pos_lose[key])
        for ele in lose_ranges:
            out_lose.write("\t".join([key, str(ele[0]), str(ele[1])]))
            out_lose.write("\n")
    out_lose.close()
    

def main():
    args = GetArgs()
    depth, cutoff, out = args.depth, args.thre, args.out
    FindRanges(depth, cutoff, out)


if __name__ == '__main__':
    main()

