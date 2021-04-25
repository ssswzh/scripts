#!/usr/bin/env python3
# Author: Zhang Siwen, Han Yue
# History:
#     20210414, first version
#     20210423, change true coordinate restore

import argparse
import pysam
import numpy as np
from pyfaidx import Fasta
import sys
from Bio.Seq import Seq
from keras.models import load_model
from pkg_resources import resource_filename
from spliceai import utils

# spliceai models
context = 10000
paths = ('models/spliceai{}.h5'.format(x) for x in range(1, 6))
models = [load_model(resource_filename('spliceai', x)) for x in paths]



def GetArgs():
    parser = argparse.ArgumentParser(description='Use SpliceAI model to predict splice site for SVs (only deletions and duplications).')
    parser.add_argument('--vcf', dest='vcf', help='SV vcf file', required=True) 
    parser.add_argument('--ref', dest='reference', help='reference fasta', required=True) 
    parser.add_argument('--version', dest='ref_version', default="grch37", help='reference version, choose from "grch37" and "grch38", default grch37')
    parser.add_argument('--flank', dest='flank', default=1000, type=int, help='flank length for deletions and duplications, default 1000, max 4999')
    parser.add_argument('--out', dest='output', help='output vcf file name', required=True) 
    args = parser.parse_args()
    return args



def ExtractSeqBySVType(reference, svtype, chrom, start, end, flanking=1000):
    fa = Fasta(reference, rebuild=False)
    if svtype == "DEL": 
        # contains neigher start or end base
        upflank = fa[chrom][start-flanking-1 : start-1].seq
        downflank = fa[chrom][end : end+flanking].seq
        #upflank = fa[chrom][start-flanking : start].seq
        #downflank = fa[chrom][end+1 : end+flanking+1].seq
    elif svtype == "DUP": 
        # contains both start and end base
        #upflank = fa[chrom][end-flanking : end].seq
        #downflank = fa[chrom][start-1 : start+flanking-1].seq
        upflank = fa[chrom][end-flanking+1 : end+1].seq
        downflank = fa[chrom][start : start+flanking].seq
    seq = upflank + downflank
    return seq


def CalculateSeqProb(seq):
    seq_matrix = utils.one_hot_encode('N'*(context//2) + seq + 'N'*(context//2))[None, :]
    seq_predict = np.mean([models[m].predict(seq_matrix) for m in range(5)], axis=0)
    '''
    acceptor_prob_list = seq_predict[0, :, 1]
    donor_prob = seq_predict[0, :, 2]
    '''
    acceptor_prob = np.amax(seq_predict[0], axis=0)[1]
    donor_prob = np.amax(seq_predict[0], axis=0)[2]
    acceptor_pos = np.where(seq_predict[0] == acceptor_prob)[0][0]
    donor_pos = np.where(seq_predict[0] == donor_prob)[0][0]
    return acceptor_prob, donor_prob, acceptor_pos, donor_pos


def GetSeqProbabilities(sequence, reference, ref_version, svtype, chrom, pos, end, flank):
    if ref_version == 'grch37':
        annotations = utils.__file__.replace("utils.py", "") + "annotations/grch37.txt"
    elif ref_version == 'grch38':
        annotations = utils.__file__.replace("utils.py", "") + "annotations/grch37.txt"
    # probabilities
    probs = []
    ann = utils.Annotator(reference, annotations)
    # format chromosome prefix
    chrom = utils.normalise_chrom(chrom, list(ann.ref_fasta.keys())[0])
    genes, strands, idxs = ann.get_name_and_strand(chrom, pos)
    if len(idxs) == 0:
        return probs
    for i in range(len(idxs)):
        # if found multiple genes, append all
        if strands[i] == "-":
            sequence = str(Seq(sequence).reverse_complement())
        # probabilities and relative positions
        acceptor_prob, donor_prob, acceptor_pos, donor_pos = CalculateSeqProb(sequence)
        # output combined sequence, upper case for pseudoexon
        if donor_pos > acceptor_pos:
            newSeq = sequence[0:acceptor_pos].lower() + sequence[acceptor_pos:donor_pos+1].upper() + sequence[donor_pos+1:].lower()
        else:
            newSeq = "--"
        # use middle position to restore relative position
        mid_pos = len(sequence)/2
        # restore true coordinate
        if svtype == "DEL":
            if strands == "+":
                if acceptor_pos < mid_pos:
                    acceptor_pos = pos-1 - flank + acceptor_pos + 1
                else:
                    acceptor_pos = end+1 - flank + acceptor_pos
                if donor_pos < mid_pos:
                    donor_pos = pos-1 - flank + donor_pos + 1
                else:
                    donor_pos = end+1 - flank + donor_pos
            elif strands == "-":
                if acceptor_pos < mid_pos:
                    acceptor_pos = end+1 + flank - acceptor_pos - 1
                else:
                    acceptor_pos = pos-1 - acceptor_pos + flank
                if donor_pos < mid_pos:
                    donor_pos = end+1 + flank - donor_pos - 1
                else:
                    donor_pos = pos-1 - donor_pos + flank
        elif svtype == "DUP":
            if strands == "+":
                if acceptor_pos < mid_pos:
                    acceptor_pos = end - flank + acceptor_pos + 1
                else:
                    acceptor_pos = pos - flank + acceptor_pos
                if donor_pos < mid_pos:
                    donor_pos = end - flank + donor_pos + 1
                else:
                    donor_pos = pos - flank + donor_pos
            elif strands == "-":
                if acceptor_pos < mid_pos:
                    acceptor_pos = pos + flank - acceptor_pos - 1
                else:
                    acceptor_pos = end + flank - acceptor_pos
                if donor_pos < mid_pos:
                    donor_pos = pos + flank - donor_pos - 1
                else:
                    donor_pos = end + flank - donor_pos
        probs.append("{}|{}|{}|{}|{}|{}|{}|{}".format(svtype, genes[i], strands[i], round(acceptor_prob,3), round(donor_prob,3), acceptor_pos, donor_pos, newSeq))
    return probs



def main():
    args = GetArgs()
    '''
    reference = "/sfs-grand-med-research/database/human/humanG1Kv37/human_g1k_v37.fasta"
    ref_version = "grch37"
    vcffile = "BM21A0457_1.chrX.vcf"
    out = "test_output.vcf"
    '''
    try:
        vcf = pysam.VariantFile(args.vcf)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant '
                    'annotation. These include acceptor probabilities (AProb), donor probabilities'
                    ' (DProb), acceptor position (APos), donor position (DPos) and sequence (Seq). '
                    'Format: SVTYPE|SYMBOL|STRAND|AProb|DProb|APos|DPos|Seq">')
    try:
        output = pysam.VariantFile(args.output, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    # read vcf
    for record in vcf:
        # deletions or duplications
        if record.info['SVTYPE'] == "DEL" or record.info['SVTYPE'] == "DUP":
            sequence = ExtractSeqBySVType(args.reference, record.info['SVTYPE'], record.chrom, record.pos, record.stop, args.flank)
            '''
            record.samples: [Call(sample=BM21A0457_1.bam, CallData(GT=0/0, DR=814, DV=6))]
            record.FORMAT: 'GT:DR:DV'
            '''
            probs = GetSeqProbabilities(sequence, args.reference, args.ref_version, record.info['SVTYPE'], record.chrom, record.pos, record.stop, args.flank)
            if len(probs) > 0:
                record.info['SpliceAI'] = probs
            output.write(record)
        else:
            output.write(record)
    output.close()



if __name__ == '__main__':
    main()


'''
How to get true genome coordinate by sequence index (pyi), flank and variant position.

DEL

1. DEL(+):21-49
pos=21, end=49, flank=10
ref: 11 12 13 14 15 16 17 18 19 20 | 50 51 52 53 54 55 56 57 58 59
seq: () () () () () () () () () () | () () () () () () () () () ()
pyi: 0  1  2  3  4  5  6  7  8  9  | 10 11 12 13 14 15 16 17 18 19
              |                                   |
              V                                   V
              14=pos-1-flank+pyi+1                54=end+1-flank+pyi

2. DEL(-):21-49
pos=21, end=49, flank=10
ref: 59 58 57 56 55 54 53 52 51 50 | 20 19 18 17 16 15 14 13 12 11
seq: () () () () () () () () () () | () () () () () () () () () ()
pyi: 0  1  2  3  4  5  6  7  8  9  | 10 11 12 13 14 15 16 17 18 19
              |                                   |
              V                                   V
              56=end+1+flank-pyi-1                16=pos-1+flank-pyi

DUP

1. DUP(+):20-50
pos=20, end=50, flank=10
ref: 41 42 43 44 45 46 47 48 49 50 | 20 21 22 23 24 25 26 27 28 29
seq: () () () () () () () () () () | () () () () () () () () () ()
pyi: 0  1  2  3  4  5  6  7  8  9  | 10 11 12 13 14 15 16 17 18 19
              |                                   |
              V                                   V
              44=end-flank+pyi+1                  24=pos-flank+pyi

2. DUP(-):20-50
pos=20, end=50, flank=10
ref: 29 28 27 26 25 24 23 22 21 20 | 50 49 48 47 46 45 44 43 42 41
seq: () () () () () () () () () () | () () () () () () () () () ()
pyi: 0  1  2  3  4  5  6  7  8  9  | 10 11 12 13 14 15 16 17 18 19
              |                                   |
              V                                   V
              26=pos+flank-pyi-1                  24=end+flank-pyi

'''


