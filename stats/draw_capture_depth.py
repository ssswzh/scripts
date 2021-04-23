#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen@grandomics.com, May 31 2019
# History:
#   July 22 2019, change a single region to a bed file and draw plot for each region
#   Feb 7 2020, add a column for 4th-end column in bed file
#   Apr 23 2020, add mean median max min value to a new outfile

import os
import argparse
import pandas as pd
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

def GetArgs():
    parser = argparse.ArgumentParser(description='Script for calculating sliding window average depth in region and draw plot.')
    parser.add_argument('--indepth', dest='indepth', help='sample.bam.depth', required=True) 
    parser.add_argument('--bed', dest='bed', help='bed file of capture region, 0-based', required=True)
    parser.add_argument('--outfile', dest='outfile', help='out file name', required=True) 
    parser.add_argument('--size', dest='size', default=2000, help='sliding window size, default=2000') 
    parser.add_argument('--step', dest='step', default=1600, help='sliding window step, default=1600') 
    args = parser.parse_args()
    return args

def ReadBedRegion(bed): # function of reading bed file, added July 22 2019
    bed = open(bed)
    bed_region = []
    for line in bed:
        ele = line.strip().split("\t")
        if len(ele) > 3: # function of bed file column number, Feb 7 2020
            info = '_'.join(ele[3:])
            bed_region.append(ele[0]+"#"+ele[1]+"#"+ele[2]+"#"+info)
        else:
            bed_region.append(ele[0]+"#"+ele[1]+"#"+ele[2]+"#")
        
    return bed_region

def ReadAndProcess(indepth, bed_region, window_size, step, outfile):
    df = pd.read_table(indepth, sep="\t", header=None, low_memory=False)
    bed_region = bed_region
    size = window_size
    step = step
    result=open(outfile, "w")
    result.write("chrom\tposition\tstart\tend\ttotal\taverage\tinfo\n")
    region_hash = {}
    for region in bed_region: # for loop to process each region, added July 22 2019
        chrom = str(region.split("#")[0])
        start = int(region.split("#")[1])
        end = int(region.split("#")[2])
        info = str(region.split("#")[3])
        df.iloc[:,0] = df.iloc[:,0].astype(str) # mandatory type conversion to string, added Jan 3 2019
        depth = df.loc[(df.iloc[:,0] == chrom) & (df.iloc[:,1] >= start) & (df.iloc[:,1] <= end)]
        # default dict value type to avoid error when use hash[key].append, added July 22 2019
        region_hash[region] = [] 
        for window_start in range(start, end, step):
            window_end = window_start+size if window_start+size < end else end
            sub = depth.loc[(depth.iloc[:,0] == chrom) & (depth.iloc[:,1] >= window_start) & (depth.iloc[:,1] < window_end)]
            total = sum(sub.iloc[:,2])
            avg = total/(window_end-window_start)
            pos = int(window_start+(window_end-window_start)/2)
            #depth.drop(depth[depth.iloc[:,1]<window_start].index, inplace=True)
            result.write(chrom + "\t" + str(pos) + "\t" + str(window_start) + "\t" + str(window_end) + "\t" + str(total) + "\t" + str(avg) + "\t" + info + "\n")
            # hash to store region as key and depth calculation as value, added July 22 2019
            region_hash[region].append([chrom, pos, window_start, window_end, total, avg, info]) 
    result.close()
    return region_hash

def DrawDepthPlot(outfile, region_hash, plot_dir):
    prefix = os.path.basename(outfile) 
    chara = "/|{}[]()>#+!$`"
    result=open(outfile+".depth_stat", "w") # count for depth stat and output in a new file
    result.write("chrom\tstart\tend\tinfo\tminDepth\tmeanDepth\tmedianDepth\tmaxDepth\n") # added Apr 23 2020
    for key in region_hash: # for loop to process each region, added July 22 2019
        key_trans = key
        for c in chara:
            key_trans = key_trans.replace(c,'.')
        outplot = plot_dir + "/" + prefix + "." + key_trans + ".pdf" # output plot dir and name, modified July 22 2019
        print(outplot)
        df = pd.DataFrame(region_hash[key], columns = ["chrom", "position", "start", "end", "total", "average", "info"])
        chrom = str(key.split("#")[0])
        start = int(key.split("#")[1])
        end = int(key.split("#")[2]) 
        info = df['info'].unique()[0] # added Apr 23 2020
        plt.figure(dpi=128,figsize=(18,6))
        plt.grid(True)
        y_max = max(df['average'])
        y_mean = df['average'].mean()
        y_median = df['average'].median()
        y_min = min(df['average'])
        
        plt.axis([start-(end-start)/20, end+(end-start)/20, 0, y_max*1.2])
        plt.vlines(int(start), 0, y_max, colors = "r", linestyles = "dashed")
        plt.vlines(int(end), 0, y_max, colors = "r", linestyles = "dashed")
        plt.hlines(y_max, start-(end-start)/20, end+(end-start)/20, colors = "blue", linestyles = "dashed")
        plt.text(end+(end-start)/20, y_max*1.03, "Max="+str(round(y_max,2)), color = 'blue', verticalalignment='center', horizontalalignment='center')
        plt.hlines(y_mean, start-(end-start)/20, end+(end-start)/20, colors = "g", linestyles = "dashed")
        plt.text(end+(end-start)/20, y_mean*1.03, "Mean="+str(round(y_mean,2)), color = 'green', verticalalignment='center', horizontalalignment='center')
        plt.hlines(y_median, start-(end-start)/20, end+(end-start)/20, colors = "r", linestyles = "dashed")
        plt.text(start-(end-start)/20, y_median*1.03, "Median="+str(round(y_median,2)), color = "red", verticalalignment='center', horizontalalignment='center')
        plt.hlines(y_min, start-(end-start)/20, end+(end-start)/20, colors = "magenta", linestyles = "dashed")
        plt.text(start-(end-start)/20, y_min*1.03, "Min="+str(round(y_min,2)), color = 'magenta', verticalalignment='center', horizontalalignment='center')
        plt.plot(df['position'], df['average']) #, c=df['average'], cmap=plt.cm.Blues, edgecolor='none', s=40)
        plt.xlabel("Coordinate")
        plt.ylabel("Avg Depth")
        plt.title("Average depth distribution " + key)    
        plt.tick_params(axis='both',which='major',labelsize=14)
        plt.savefig(outplot, bbox_inches='tight')

        result.write('\t'.join([chrom,str(start),str(end),info,str(int(y_min)),str(round(y_mean,2)),str(round(y_median,2)),str(int(y_max))]) +"\n")
    result.close()

def main():
    args = GetArgs()
    infile = args.indepth
    outfile = args.outfile
    bed = args.bed #X:31115345-33377726
    window_size = int(args.size)
    step = int(args.step)
    # put all plots inside the folder as outfile, added July 22 2019
    plot_dir = outfile+"_plot/"
    print(plot_dir)
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)
    bed_region = ReadBedRegion(bed) # added July 22 2019
    region_hash = ReadAndProcess(infile, bed_region, window_size, step, outfile)
    DrawDepthPlot(outfile, region_hash, plot_dir)

if __name__=="__main__":
    main()
