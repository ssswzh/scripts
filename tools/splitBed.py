#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2021/12/29
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2021/12/29
# @ChangeLog
#     20211229, first version

import argparse

def getArgs():
    parser = argparse.ArgumentParser(description="chop bed file into customized size and step")
    parser.add_argument('--bed', help='input bed file', action='store', dest='bed', required=True)
    parser.add_argument('--size', help='region size in bed file', action='store', dest='size', required=True)
    parser.add_argument('--step', help='region step, default equals to region size', action='store', required=False)
    parser.add_argument('--out', help='out bed file', action='store', dest='out', required=True)
    args = parser.parse_args()
    return args

def splitBed(bedfile, size, step, outfile):
    out = open(outfile, "w")
    with open(bedfile) as bedinfo:
        for line in bedinfo:
            regions = []
            record = line.rstrip().split("\t")
            chrom, start, end = record[0], int(record[1]), int(record[2])
            if len(record) > 3:
                otherinfo = "\t".join(record[3:])
            else:
                otherinfo = ""
            if size != step:
                for i in range(1, size, step):
                    for j in range(i, end, size):
                        if j+size > end:
                            regions.append([chrom, j, end])
                        else:
                            regions.append([chrom, j, j+size])
            else:
                for i in range(1, end, size):
                    if i+size > end:
                        regions.append([chrom, i, end])
                    else:
                        regions.append([chrom, i, i+size])
            # sort by start pos
            regions = sorted(regions, key = lambda x: x[1])
            for i in regions:
                if otherinfo:
                    out.write("\t".join([i[0], str(i[1]), str(i[2]), otherinfo]) + "\n")
                else:
                    out.write("\t".join([i[0], str(i[1]), str(i[2])]) + "\n")
    out.close()

def main():
    """ """
    args = getArgs()
    if not args.step:
        args.step = args.size
    splitBed(args.bed, int(args.size), int(args.step), args.out)


if __name__ == "__main__":
    main()

