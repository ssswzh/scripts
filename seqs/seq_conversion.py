#!/usr/bin/env python3
# coding: utf-8

import os
import argparse
from Bio.Seq import Seq
from Bio.SeqUtils import GC


def get_args():
    parser = argparse.ArgumentParser(description="fast sequence transform", 
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    parser.add_argument("-f", required=True, type=str, 
                                help="Tab-delimited format file OR sequence")
    parser.add_argument("-c", required=False, type=int, default=1, 
                                help="1-based column number of sequence to be transformed")
    parser.add_argument("-m", required=False, type=str, default="rev_com", 
                                choices=["rev","com","rev_com"], 
                                help="reverse|complement|reverse_complement")
    parser.add_argument("-gc", required=False, action="store_true", 
                                help="also get sequence GC to new column")
    usage = '''Usage:
    %(prog)s --in seq.tsv --col 1 --mode [rev|com|rev_com]
    %(prog)s --in GTTTAATTGAGTTGTCATATGTTAATAACGGTAT
    rev: reverse
    com: complement
    rev_com: reverse_complement
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


def seq_transform(sequence, mode="rev_com"):
    if mode == "rev_com":
        return str(Seq(sequence).reverse_complement())
    if mode == "rev":
        return str(sequence[::-1])
    if mode == "com":
        return str(Seq.complement(Seq(sequence)))
    else:
        return None


def main():
    args = get_args()
    if os.path.isfile(str(args.f)):
        infile = open(args.f)
        column = int(args.c) - 1
        for line in infile:
            seq = seq_transform(line.strip().split("\t")[column], mode=args.m)
            if args.gc:
                print(line.strip() + "\t" + seq + "\t" + str(GC(seq)))
            else:
                print(line.strip() + "\t" + seq)
    else:
        seq = str(args.f)
        if args.gc:
            print(seq_transform(seq, mode=args.m) + "\t" + str(GC(seq)))
        else:
            print(seq_transform(seq, mode=args.m))
    return None


if __name__ == '__main__':
    main()
