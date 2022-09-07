#!/usr/bin/env python3
import Levenshtein
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description = "script for Levenshtein distance and ratio calculation for fasta sequences")
parser.add_argument('-i', help='input fasta file for Levenshtein distance calculation', action='store', dest='input', required=True)
parser.add_argument('-o', help='output prefix', action='store', dest='output', required=True)
args = parser.parse_args()

#fastafile = "highrisk.E7.fasta"
seqs = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
'''
seqs['NC_001526.4:7603-7900']
SeqRecord(seq=Seq('ATGCATGGAGATACACCTACATTGCATGAATATATGTTAGATTTGCAACCAGAG...TAA'), id='NC_001526.4:7603-7900', name='NC_001526.4:7603-7900', description='NC_001526.4:7603-7900', dbxrefs=[])
'''

out = open(args.output, "w")
keys = list(seqs.keys())
out.write("Ratio|Distance\t" + "\t".join(keys) + "\n")
for i in range(0,len(keys)):
    matrix = []
    for j in range(0,len(keys)):
        r = str(Levenshtein.seqratio(str(seqs[keys[i]].seq), str(seqs[keys[j]].seq)))
        d = str(Levenshtein.distance(str(seqs[keys[i]].seq), str(seqs[keys[j]].seq)))
        if j < i:
            matrix.append(r)
        if j == i:
            matrix.append("1")
        if j > i:
            matrix.append(d)
    out.write(keys[i] + "\t" + "\t".join(matrix) + "\n")
out.close()

