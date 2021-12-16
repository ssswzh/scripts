#!/usr/bin/env python3
# coding: utf-8

from Bio.Seq import Seq
from Bio.SeqUtils import GC
import sys

if len(sys.argv) == 1 or len(sys.argv) !=3:
    print("\npython3 {} <file> <0based_column_to_be_reversed_complement> > out.file\n".format(sys.argv[0]))
    sys.exit()
elif sys.argv[1] == "-h" or  sys.argv[1] == "--help":
    print("\npython3 {} <file> <0based_column_to_be_reversed_complement> > out.file\n".format(sys.argv[0]))
    sys.exit()

infile = open(sys.argv[1])
column = int(sys.argv[2])
#outfile = open(sys.argv[1]+".revcomp", "w")

for line in infile:
    ele = line.strip().split("\t")
    seq = ele[column]
    revcomp = str(Seq(seq).reverse_complement())
    #outfile.write(line.strip()+"\t"+revcomp+"\n")
    print(line.strip()+"\t"+revcomp+"\t"+str(GC(revcomp)))

infile.close()
#outfile.close()
#print("Wrote to {} \n".format(sys.argv[1]+".revcomp"))
