#!/usr/bin/env python3
# coding: utf-8
# Author:	jinhongshuai
# Date & Time:	2019-06-24

import sys,re,os,pysam,gzip
import argparse
import pandas as pd

def MKDIR(*dirs):
	for dir in dirs:
		if os.path.exists(dir):
			pass
		else:
			os.mkdir(dir)

class FinalStat(object):

	def __init__(self):
		self.capture_bed = None
		self.outdir = None
		self.prefix = None

	def Getopt(self):
		parser = argparse.ArgumentParser(description='Bam depth stat')
		parser.add_argument('--prefix', dest='prefix', help='prefix,str')
		parser.add_argument('--MapStat', dest='MapStat', help='MapStat,file')
		parser.add_argument('--bam', dest='bam', help='bam,file')
		parser.add_argument('--capture_bed', dest='capture_bed', help='bed,file')
		parser.add_argument('--mosdepth_regions', dest='mosdepth_regions', help='regions,array')
		parser.add_argument('--mosdepth_thresholds', dest='mosdepth_thresholds', help='thresholds,array')
		parser.add_argument('--mosdepth_per_base', dest='mosdepth_per_base', help='per-base,array')
		parser.add_argument('--outdir', dest='od', help='output,dir')
		opt = parser.parse_args()
		
		self.MapStat=opt.MapStat
		if opt.capture_bed:
			self.capture_bed=opt.capture_bed
		self.bam=opt.bam
		self.mosdepth_regions_list=(opt.mosdepth_regions).split(",")
		self.mosdepth_thresholds_list=(opt.mosdepth_thresholds).split(",")
		self.mosdepth_per_base_list=(opt.mosdepth_per_base).split(",")
		self.outdir=opt.od
		MKDIR(self.outdir)
		self.prefix=opt.prefix
		
	def ReadMapStat(self):
		dictMapinfo={}
		for line in open(self.MapStat):
			if line.startswith("#"):
				continue
			arr=line.strip().split('\t')
			if arr[0]=='SN':
				if arr[1]=='sequences:':
					total_reads=arr[2]
					dictMapinfo['total_reads']=total_reads
				if arr[1]=='total length:':
					total_base=arr[2]
					dictMapinfo['total_base']=total_base
				if arr[1]=='reads mapped:':
					mapped_reads=arr[2]
					dictMapinfo['mapped_reads']=mapped_reads
				if arr[1]=='bases mapped:':
					mapped_base=arr[2]
					dictMapinfo['mapped_base']=mapped_base
		return dictMapinfo
		
	def CaptureStat(self):
		dictMappedReads={}
		
		df=pd.read_csv((self.mosdepth_per_base_list)[-1],compression ='gzip',sep="\t",header=None,error_bad_lines = False,low_memory=False)
		sam=pysam.AlignmentFile(self.bam,'rb')
		
		for line in open(self.capture_bed):
			dictReads={}
			CHR,START,END=line.strip().split('\t')[:3]
			location=CHR+'_'+START+'_'+END
			for r in sam.fetch('%s'%CHR,int(START),int(END)):
				dictReads[r.query_name]=''
			MappedReads=len(dictReads)
			
			loca_df=df[(df[0]==CHR)&(df[2]>=int(START))&(df[1]<=int(END))]
			MappedBase=((loca_df[2]-loca_df[1])*loca_df[3]).sum()
			
			dictMappedReads[location]=[MappedReads,MappedBase]
			
		sam.close()
		return dictMappedReads
	
	def OutStat(self,dictMapinfo,dictMappedReads):
		f=open(self.outdir+'/'+self.prefix+'.stat.xls','w')
		f.write('#total_reads	total_base	mapped_reads	mapped_reads_ratio	mapped_base	mapped_base_ratio\n')
		f.write(str(dictMapinfo['total_reads'])+'\t'+str(dictMapinfo['total_base'])+'\t'+str(dictMapinfo['mapped_reads'])+'\t'+str(float(dictMapinfo['mapped_reads'])/float(dictMapinfo['total_reads'])*100)+'\t'+str(dictMapinfo['mapped_base'])+'\t'+str(float(dictMapinfo['mapped_base'])/float(dictMapinfo['total_base'])*100)+'\n')
		f.write('#position	AverageDepth	MappedReads	MappedReads_ratio	MappedBase	MappedBase_ratio\n')
		dictAverageDepth={}
		for line in gzip.open(self.mosdepth_regions_list[-1]):
			line=line.decode()##python3
			arr=line.strip().split("\t")
			CHR,S,E=arr[:3]
			DEPTH=arr[-1]
			dictAverageDepth[CHR+'_'+S+'_'+E]=DEPTH
		for pos in dictMappedReads:
			f.write(pos+'\t'+dictAverageDepth[pos]+"\t"+str(dictMappedReads[pos][0])+'\t'+str(float(dictMappedReads[pos][0])/float(dictMapinfo['total_reads'])*100)+'\t'+str(dictMappedReads[pos][1])+'\t'+str(float(dictMappedReads[pos][1])/float(dictMapinfo['total_base'])*100)+'\n')
		
		f.close()
		
	def run(self):
		self.Getopt()
		dictMappedReads=self.CaptureStat()
		if self.capture_bed:
			dictMapinfo=self.ReadMapStat()
			self.OutStat(dictMapinfo,dictMappedReads)
			
if __name__ == "__main__":
	runFinalStat=FinalStat()
	runFinalStat.run()
