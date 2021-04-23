#!/usr/bin/env python3

import argparse
import gzip
import random
from readfq import readfq


FileTypeFA = (".fasta", ".fa", ".FASTA", ".FA")
FileTypeFAGZ = (".fasta.gz", ".fa.gz", ".FASTA.GZ", ".FA.GZ", ".FASTA.gz", ".FA.gz", ".fasta.GZ", ".fa.GZ")
FileTypeFQ = (".fastq", ".fq", ".FASTQ", ".FQ")
FileTypeFQGZ = (".fastq.gz", ".fq.gz", ".FASTQ.GZ", ".FQ.GZ", ".FASTQ.gz", ".FQ.gz", ".fastq.GZ", ".fq.GZ")

def GetArgs():
    parser = argparse.ArgumentParser(description='Split FASTA/FASTQ file each read by random length, not split file by random number of reads.')
    parser.add_argument('--input', dest='input', help='input FASTA/FASTQ file', required=True) 
    parser.add_argument('--min', dest='minlength', default=2500, type=int, help='the lower threshold for random split length, default 2500')
    parser.add_argument('--max', dest='maxlength', default=3500, type=int, help='the upper threshold for random split length, default 3500')
    parser.add_argument('--residuel', dest='residuel', default=500, type=int, help='the cutoff to keep residuel sequence of splitted reads, default 500') 
    parser.add_argument('--out', dest='output', help='output file name', required=True) 
    args = parser.parse_args()
    return args


def CheckFileFileType(input):
    if input.endswith(FileTypeFA):
        FileType, gz = "FA", False
    elif input.endswith(FileTypeFAGZ):
        FileType, gz = "FA", True
    elif input.endswith(FileTypeFQ):
        FileType, gz = "FQ", False
    elif input.endswith(FileTypeFQGZ):
        FileType, gz = "FQ", True
    else:
        print("Please provide file ends with:")
        print(FileTypeFA, FileTypeFAGZ)
        print(FileTypeFQ, FileTypeFQGZ)
    return FileType, gz


def RandomLengthGeneration(seq, minlength, maxlength, residuel):
    position = [0]
    while True:
        if position[-1] < len(seq):
            l = random.randint(minlength, maxlength)
            if l + position[-1] <= len(seq) and l >= residuel:
                position.append(position[-1] + l)
                continue
            elif l + position[-1] >= len(seq) and len(seq) - position[-1] >= residuel:
                position.append(len(seq))
                break
            else:
                break
        else:
            break
    return position


# split sequences by two neighbouring coordinates
def SplitBYLength(seq, position):
    sequences = [seq[position[i]:position[i+1]] for i in range(0, len(position)-1)]
    return sequences


def SplitFastqSequences(infile, minlength, maxlength, residuel, outfile, FileType):
    position = []
    read_name = ""
    seq = ""
    split_seq = []
    split_qual = []
    for read_name, sequence, quality in readfq(infile):
        position = RandomLengthGeneration(sequence, minlength, maxlength, residuel)
        split_seq = SplitBYLength(sequence, position)
        if FileType == "FA":
            for i in range(0,len(split_seq)):
                outfile.write(">"+read_name+"/"+str(i+1)+"\n")
                outfile.write(split_seq[i]+"\n")
        if FileType == "FQ":
            split_qual = SplitBYLength(quality, position)
            for i in range(0,len(split_qual)):
                outfile.write("@"+read_name+"/"+str(i+1)+"\n")
                outfile.write(split_seq[i]+"\n")
                outfile.write("+\n")
                outfile.write(split_qual[i]+"\n")


def main():
    args = GetArgs()
    FileType, gz = CheckFileFileType(args.input)
    if gz:
        with gzip.open(args.input, 'rb') as infile, gzip.open(args.output, 'wb') as outfile:
            SplitFastqSequences(infile, args.minlength, args.maxlength, args.residuel, outfile, FileType)
    else:
        with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
            SplitFastqSequences(infile, args.minlength, args.maxlength, args.residuel, outfile, FileType)
    infile.close()
    outfile.close()
    
if __name__ == '__main__':
    main()
