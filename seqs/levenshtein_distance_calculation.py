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

outratio = open(args.output + ".ratio", "w")
outdist = open(args.output + ".distance", "w")
keys = list(seqs.keys())
outratio.write("Ratio\t" + "\t".join(keys) + "\n")
outdist.write("Distance\t" + "\t".join(keys) + "\n")
for i in range(0,len(keys)):
    ratios = []
    distances = []
    for j in range(0,len(keys)):
        ratios.append(str(Levenshtein.seqratio(str(seqs[keys[i]].seq), str(seqs[keys[j]].seq))))
        distances.append(str(Levenshtein.distance(str(seqs[keys[i]].seq), str(seqs[keys[j]].seq))))
    outratio.write(keys[i] + "\t" + "\t".join(ratios) + "\n")
    outdist.write(keys[i] + "\t" + "\t".join(distances) + "\n")
    #outratio.write(keys[i] + "\t"*i + "\t".join(ratios) + "\n")
    #outdist.write(keys[i] + "\t"*i + "\t".join(distances) + "\n")
outratio.close()
outdist.close()

