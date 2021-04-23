#!/export/software/conda/miniconda3/bin/python
# coding: utf-8
# Author: zhangsiwen@grandomics.com, Sept 5th 2019
# History:
#   Sept 6th 2019, add bin segmentation function to avoid large HTML file

import pysam
import argparse
import pandas as pd
import plotly.graph_objs as go


def GetArgs():
    parser = argparse.ArgumentParser(description='Script for calculating read count and position depth in region and draw plot in HTML.')
    parser.add_argument('--bam', dest='bam', help='sample.bam', required=True) 
    parser.add_argument('--bed', dest='bed', help='bed file of capture region, chr\\tstart\\tend\\tname', required=True)
    parser.add_argument('--outprefix', dest='out', help='out file prefix', required=True) 
    parser.add_argument('--overlap', dest='overlap', default=0.5, help='sliding window size, default=0.5')
    parser.add_argument('--exon', dest='exon', default=10, help='bin size of exon, default=10')
    parser.add_argument('--intron', dest='intron', default=1000, help='bin size of intron, default=1000')
    parser.add_argument('--whole', dest = 'whole', action = 'store_true', help = 'use this if want depth in every position in bed')
    args = parser.parse_args()
    return args


def ReadBedRegion(bed):
    bed_file = open(bed)
    region = []
    # region from bed
    print("Start read bed file.") 
    for line in bed_file:
        ele = line.strip().split("\t")
        region.append([ele[0],ele[1],ele[2],ele[3]])
    bed_file.close()
    print("Read bed file and store region --> Done.")
    return region


def RegionReadCount(bam, region, overlap_ratio, out_prefix):
    bam_file = pysam.AlignmentFile(bam, "rb" )
    final_count = {}
    overlap_ratio = int(overlap_ratio)
    out_readcount = open(out_prefix+".read_count","w")
    out_readcount.write("Chrom\tStart\tEnd\tName\tLength\tReadCount\n")
    # read count part
    print("Start calculating read count in each region.")
    for bed in region:
        chrom = bed[0]
        start = int(bed[1])
        end = int(bed[2])
        name = bed[3]
        length = end-start+1
        min_overlap = overlap_ratio*length
        read_count = 0
        for aln in bam_file.fetch(chrom, start, end):
            # calculate read count if overlap > ratio*length
            overlap = min(aln.reference_end, end) - max(aln.reference_start, start)
            if overlap >= min_overlap:
                read_count +=1
        final_count[name] = read_count
        out_readcount.write(chrom+"\t"+str(start)+"\t"+str(end)+"\t"+name+"\t"+str(length)+"\t"+str(read_count)+"\n")
    bam_file.close()
    out_readcount.close()
    print("Read count calculation --> Done.")
    #return final_count


def PositionDepthDataframe(bam, region, out_prefix):
    out_depth = open(out_prefix+".depth","w")
    out_depth.write("Chrom\tStart\tEnd\tName\tLength\tDepth\n")
    # position depth plot part
    bamfile = pysam.AlignmentFile(bam, "rb" )
    print("Start parsing bam file for each position in bed.")
    exon_x = []
    exon_y = []
    exon_text = []
    intron_x = []
    intron_y = []
    intron_text = []
    for bed in region:
        chrom = bed[0]
        start = int(bed[1])
        end = int(bed[2])
        name = bed[3]
        length = end-start+1
        if "exon" in name:
            for pileupcolumn in bamfile.pileup(chrom, start, end):
                pos = pileupcolumn.pos
                depth = pileupcolumn.nsegments
                #depth = sum(bam_file.count_coverage(chrom, pos, pos+1, quality_threshold = 0)[i][0] for i in range(0,4))
                if pos >= start and pos <= end:
                    exon_x.append(pos)
                    exon_y.append(depth)
                    exon_text.append(name)
                    out_depth.write(chrom+"\t"+str(pos)+"\t"+str(pos)+"\t"+name+"\t"+str(length)+"\t"+str(depth)+"\n")
        elif "intron" in name:
            for pileupcolumn in bamfile.pileup(chrom, start, end):
                pos = pileupcolumn.pos
                depth = pileupcolumn.nsegments
                #depth = sum(bam_file.count_coverage(chrom, pos, pos+1, quality_threshold = 0)[i][0] for i in range(0,4))
                if pos >= start and pos <= end:
                    intron_x.append(pos)
                    intron_y.append(depth)
                    intron_text.append(name)
                    out_depth.write(chrom+"\t"+str(pos)+"\t"+str(pos)+"\t"+name+"\t"+str(length)+"\t"+str(depth)+"\n")
    print("Bam file parsing --> Done.")
    bamfile.close()
    out_depth.close()
    exon_df = pd.DataFrame({'exon_x':exon_x, 'exon_y':exon_y, 'exon_text':exon_text})
    intron_df = pd.DataFrame({'intron_x':intron_x, 'intron_y':intron_y, 'intron_text':intron_text})
    return exon_df, intron_df


def BinDepthDataframe(bam, region, out_prefix, exon_size, intron_size):
    out_suffix = ".exon"+str(exon_size)+"_intron"+str(intron_size)
    out_depth = open(out_prefix+out_suffix+".depth", "w")
    out_depth.write("Chrom\tPosition\tStart-End\tName\tLength\tDepth\n")
    bamfile = pysam.AlignmentFile(bam, "rb" )
    print("Start parsing bam file for each region in bed.")
    print("Bin size: "+out_suffix)
    exon_x = []
    exon_y = []
    exon_text = []
    intron_x = []
    intron_y = []
    intron_text = []
    for bed in region:
        chrom = bed[0]
        start = int(bed[1])
        end = int(bed[2])
        name = bed[3]
        length = end-start+1
        x = []
        y = []
        text = []
        bin_size = exon_size if "exon" in name else intron_size
        for bin_start in range(start, end, bin_size):
            bin_end = bin_start + bin_size if bin_start + bin_size < end else end
            for pileupcolumn in bamfile.pileup(chrom, bin_start, bin_end):
                pos = pileupcolumn.pos
                depth = pileupcolumn.nsegments
                #depth = sum(bam_file.count_coverage(chrom, pos, pos+1, quality_threshold = 0)[i][0] for i in range(0,4))
                if pos >= bin_start and pos <= bin_end:
                    y.append(depth)
            bin_x_mean = round((bin_end+bin_start)/2,1)
            bin_y_mean = round(sum(y)/len(y))
            if "exon" in name:
                exon_x.append(bin_x_mean)
                exon_y.append(bin_y_mean)
                exon_text.append(name)
                out_depth.write(chrom+"\t"+str(bin_x_mean)+"\t"+str(bin_start)+"-"+str(bin_end)+"\t"+name+"\t"+str(bin_end-bin_start)+"\t"+str(bin_y_mean)+"\n")
            elif "intron" in name:
                intron_x.append(bin_x_mean)
                intron_y.append(bin_y_mean)
                intron_text.append(name)
                out_depth.write(chrom+"\t"+str(bin_x_mean)+"\t"+str(bin_start)+"-"+str(bin_end)+"\t"+name+"\t"+str(bin_end-bin_start)+"\t"+str(bin_y_mean)+"\n")
    print("Bam file parsing --> Done.")
    bamfile.close()
    out_depth.close()
    exon_df = pd.DataFrame({'exon_x':exon_x, 'exon_y':exon_y, 'exon_text':exon_text})
    intron_df = pd.DataFrame({'intron_x':intron_x, 'intron_y':intron_y, 'intron_text':intron_text})
    return exon_df, intron_df, out_suffix


def DepthToHTML(exon_df, intron_df, out_prefix):
    fig = go.Figure()
    for key in exon_df['exon_text'].unique():
        exon_x = list(exon_df[exon_df['exon_text'].str.match(key)]['exon_x'])
        exon_y = list(exon_df[exon_df['exon_text'].str.match(key)]['exon_y'])
        exon_text = list(exon_df[exon_df['exon_text'].str.match(key)]['exon_text'])
        fig.add_trace(go.Scatter(x=exon_x, y=exon_y, text=exon_text, mode = "lines+markers", name = key))
    for key in intron_df['intron_text'].unique():
        intron_x = list(intron_df[intron_df['intron_text'].str.match(key)]['intron_x'])
        intron_y = list(intron_df[intron_df['intron_text'].str.match(key)]['intron_y'])
        intron_text = list(intron_df[intron_df['intron_text'].str.match(key)]['intron_text'])
        fig.add_trace(go.Scatter(x=intron_x, y=intron_y, text=intron_text, mode = "lines+markers", name = key))
    fig.update_xaxes(title_text='Coordinate', title_font=dict(size=20), tickfont=dict(size=16), showgrid=True)
    fig.update_yaxes(title_text='Depth', title_font=dict(size=20), tickfont=dict(size=16), showgrid=True)
    out_html = open(out_prefix+".depth.html","w")
    out_html.write(fig.to_html())
    out_html.close()
    print("HTML writting --> Done.")


def main():
    args = GetArgs()
    region = ReadBedRegion(args.bed) #ReadBedRegion("DMD_capture.bed.head20")
    RegionReadCount(args.bam, region, args.overlap, args.out) #RegionReadCount("D22.bam", region, 0.5, "D22.out")
    if args.whole:
        exon_df, intron_df = PositionDepthDataframe(args.bam, region, args.out)
        DepthToHTML(exon_df, intron_df, args.out)
    exon_df, intron_df, out_suffix = BinDepthDataframe(args.bam, region, args.out, args.exon, args.intron)
    DepthToHTML(exon_df, intron_df, args.out+out_suffix)


if __name__ == "__main__":
    main()
