#!/usr/bin/env python3
# Author: siwen
# History:
#     20210331, manuscripy
#     20210407, only compare regions inside bed file, not full-length read

import re
import argparse


def GetArgs():
    parser = argparse.ArgumentParser(description='Compare full length read map info to splited read map info, give correctly mapped read number.')
    parser.add_argument('--full', dest='full', help='full-length read map info', required=True)
    parser.add_argument('--split', dest='split', help='splited read map info', required=True)
    parser.add_argument('--bed', dest='bed', help='bed file', required=True)
    parser.add_argument('--distance', dest='distance', default=100, type=int, help='the distance between full-length read mapped position and splited read mapped position,default 100')
    parser.add_argument('--out', dest='output', help='output file name', required=True)
    args = parser.parse_args()
    return args


def BEDparsing(bedfile):
    bed = []
    while open(bedfile) as b:
        for line in b:
            bed.append(line.strip().split("\t"))
    return bed


def checkKEYinfo(read, d, info):
    # read: read_name
    # d: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    chrom, chrstart, chrend = info
    if read in d:
        if chrom in d[read]:
            d[read][chrom].append([chrstart, chrend])
        else:
            d[read][chrom] = [[chrstart, chrend]]
    else:
        d[read] = {chrom: [[chrstart, chrend]]}


def AlnPosFlatten(dic, distance):
    # dic: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    # same as d in checkKEYinfo(read, d, info)
    newdic = {}
    for read in dic:
        newdic[read] = {}
        for chrom in dic[read]:
            l = dic[read][chrom]
            # [[63074325,63077244], [63077244,63083955], [63083960,63085364], [77146724,77147929], [107126268,107128311], [107128311,107131313], [107131313,107134585]]
            l = sorted(l, key=lambda x:x[0])
            newl = []
            start, end = l[0]
            for r in l[1:]:
                if start - distance <= r[0] <= end + distance:
                    if r[1] >= end + distance:
                        end = r[1]
                else:
                    newl.append([start, end])
                    start, end = r
            newl.append([start, end]) # the last region stored in variable start and end
            newdic[read][chrom] = newl
    return newdic


def AlnInfo(infile, distance, bed):
    aln = {}
    for line in open(infile,"r"):
        if line.startswith("ReadID"):
            continue
        elem = line.strip().split("\t")
        read_name, chrom, chrstart, chrend, strand = elem[0], elem[4], int(elem[5]), int(elem[6]), elem[9]
        readstart, readend, aln_length, read_length = int(elem[1]), int(elem[2]), int(elem[3]), int(elem[7])
        #generate read_info 20210407
        for 
        read_info = [chrom, chrstart, chrend]
        if read_name.endswith("ccs"): # full length read aln
            checkKEYinfo(read_name, aln, read_info)
        else: # splited read aln
            repos = re.search(r"/[0-9]+$", read_name).span()[0]
            read_id = read_name[0:repos]
            checkKEYinfo(read_id, aln, read_info)
    # flatten list of positions
    aln = AlnPosFlatten(aln, distance)
    return aln


def AlnCompare(full_aln, split_aln, distance, outfile):
    # full_aln: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    # newl = [[63074325, 63085364], [77146724, 77147929], [107126268, 107134585]]
    # a = [[63074325, 63085354], [77146724, 77147919], [107126268, 107134785]]
    # b = [[63074325, 63085354], [107126268, 107134785]]
    out = open(outfile, "w")
    out.write("ReadID\tFull-length Read Align Info\tSplited Read Align Info\tDecision\n")
    full_key_number = len(full_aln.keys())
    split_key_number = len(split_aln.keys())
    correct = 0 # cover true chrom region
    half_correct = 0 # cover sub-region of true chrom region
    partial_correct = 0 # cover one or more of true chrom regions, also mapped to other chrom or other region
    incorrect = 0 # unmapped or none of the chrom or region correct


def main():
    args = GetArgs()
    bed = BEDparsing(args.bed)
    full_aln = AlnInfo(args.full, args.distance, bed)
    split_aln = AlnInfo(args.split, args.distance, bed)
    AlnCompare(full_aln, split_aln, args.distance, args.output)


if __name__ == '__main__':
    main()