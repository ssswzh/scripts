#!/usr/bin/python
# coding: utf-8
# Author: zhangsiwen, June 13 2019
# Reference:
#   https://www.jb51.net/article/142486.htm
# History:
#   20190801, add alignment length, start and end in read
#   20190926, add shell script to get bases between two splits of one read
#   20210331, check if read is mapped

import sys
import argparse
import pysam
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def GetArgs():
    parser = argparse.ArgumentParser(description='Obtain mapping start and end position and draw a plot.')
    parser.add_argument('--bam', dest='bam', help='sample.bam', required=True) 
    parser.add_argument('--out', dest='out', help='out file prefix', required=True) 
    parser.add_argument('--seq', dest='seq', help='choose to print sequences', default=False, action='store_true')
    parser.add_argument('--plot', dest='plot', help='choose to draw alignment position plot', default=False, action='store_true')
    args = parser.parse_args()
    return args

def ParseBamFile(bam, out, seq):
    bam_file = pysam.AlignmentFile(bam, "rb")
    out_file = open(out, "w")
    if seq:
        out_file.write("ReadID\tRead_start\tRead_end\tAlignment_length\tChr\tAlignment_start\tAlignment_end\tRead_length\tMapping_Quality\tStrand\tSequence\n") # modified Aug 1 2019
    else:
        out_file.write("ReadID\tRead_start\tRead_end\tAlignment_length\tChr\tAlignment_start\tAlignment_end\tRead_length\tMapping_Quality\tStrand\n") # modified Aug 1 2019
    for line in bam_file:
        if line.is_reverse:
            strand = "-"
        else:
            strand = "+"
        if line.is_unmapped: # added 20210331
            continue
        else:
            if seq:
                out_file.write(str(line.qname)+"\t"+str(line.query_alignment_start+1)+"\t"+str(line.query_alignment_end)+"\t"+str(line.query_alignment_length)+"\t"+str(line.reference_name)+"\t"+str(line.reference_start+1)+"\t"+str(line.reference_end)+"\t"+str(line.query_length)+"\t"+str(line.mapping_quality)+"\t"+strand+"\t"+str(line.query_sequence)+"\n") # modified Aug 1 2019
            else:
                out_file.write(str(line.qname)+"\t"+str(line.query_alignment_start+1)+"\t"+str(line.query_alignment_end)+"\t"+str(line.query_alignment_length)+"\t"+str(line.reference_name)+"\t"+str(line.reference_start+1)+"\t"+str(line.reference_end)+"\t"+str(line.query_length)+"\t"+str(line.mapping_quality)+"\t"+strand+"\n") # modified Aug 1 2019
    bam_file.close()
    out_file.close()
    

def AutoLabel(rects):
    for rect in rects:
        height = rect.get_height()
        pos = rect.get_x()
        plt.text(rect.get_x()-rect.get_width()*5, 1.03*height, '%s' % int(pos))

def DrawDepthPlot(outfile):
    df = pd.read_table(outfile, sep="\t")
    # extract position and number
    start = df['Alignment_start']
    start_x = []
    start_num = []
    for i in np.unique(start):
        start_x.append(i)
        start_num.append(sum(start==i))
    end = df['Alignment_end']
    end_x = []
    end_num = []
    for i in np.unique(end):
        end_x.append(i)
        end_num.append(sum(end==i))
    
    # start position plot
    plt.figure(dpi=128,figsize=(12,6))
    plt.grid(True)
    plt.axis([min(start_x)-(max(start_x)-min(start_x))/20, max(start_x)+(max(start_x)-min(start_x))/20, 0, max(start_num)*1.1])
    bar_width = int((max(start_x)-min(start_x))/len(start_x))/5
    start_plot = plt.bar(start_x, start_num, width=bar_width)
    AutoLabel(start_plot)
    plt.xlabel("Position")
    plt.ylabel("Read Number")
    plt.title("Alignment Start Position")    
    plt.tick_params(axis='both',which='major',labelsize=14)
    plt.savefig(outfile+".start_pos.pdf", bbox_inches='tight')
    
    # end position plot
    plt.figure(dpi=128,figsize=(12,6))
    plt.grid(True)
    plt.axis([min(end_x)-(max(end_x)-min(end_x))/20, max(end_x)+(max(end_x)-min(end_x))/20, 0, max(end_num)*1.1])
    bar_width = int((max(end_x)-min(end_x))/len(end_x))/5
    end_plot = plt.bar(end_x, end_num, width=bar_width)
    AutoLabel(end_plot)
    plt.xlabel("Position")
    plt.ylabel("Read Number")
    plt.title("Alignment End Position")    
    plt.tick_params(axis='both',which='major',labelsize=14)
    plt.savefig(outfile+".end_pos.pdf", bbox_inches='tight')


def main():
    args = GetArgs()
    ParseBamFile(args.bam, args.out, args.seq)
    if args.plot:
        DrawDepthPlot(args.out)

if __name__=="__main__":
    main()
