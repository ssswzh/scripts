#!/usr/bin/python3
# coding: utf-8

"""
SCAs 标准流程 包括PB、ONT数据
"""

__author__ = 'Jin hongshuai'
__email__ = 'jinhongshuai@sina.cn'
__version__ = '0.0.1'
__status__ = 'Dev'

import sys,os,re,time,gzip,json
import argparse
import logging
#import pandas as pd
import subprocess ##python3

def MKDIR(*dirs):
	for dir in dirs:
		if os.path.exists(dir):
			pass
		else:
			os.mkdir(dir)
			
def TIME():
	mytime=time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
	print mytime
	
def ABS(file):
	return os.path.abspath(file)
	
	
	
def get_args():
	parser = argparse.ArgumentParser(description='Bam depth stat')
	parser.add_argument('-v','--vcf', metavar='vcf', help='')
	parser.add_argument('-g','--pileup' , metavar='pileup',help='')
	parser.add_argument('-p','--prefix', metavar='prefix', help='')
	parser.add_argument('-o','--outdir' , metavar='outdir' , help='')
	opt = parser.parse_args()
	return opt
	
def read_vcf(vcf):
	dictvcf={}
	pos_list=[]
	for line in open(vcf):
		if line.startswith("#"):
			continue
		chr,pos,point,ref,alt=line.strip().split("\t")[:5]
		info=line.strip().split("\t")[9].split(":")[0]
		if not re.search("\|",info):
			continue
		type1,type2=info.split("|")
		if type1=="0":
			dictvcf[pos]=[ref,alt]
		if type1=="1":
			dictvcf[pos]=[alt,ref]
		pos_list.append(pos)
		
	return dictvcf,pos_list
	
def process_pileup(chrom,ref_coordinate,refbase,depth,mapinfo,quality,reads_coordinate,reads):
	
	mapinfo_list=[0 for n in range(int(depth))]
	templist=list(mapinfo)
	for index in range(int(depth)):
		if templist[0] in ['^']:                          ## ^]. 起始碱基 ^后面的字符ascii码减去33代表比对质量
			sign=("").join(templist[0:3])
			del(templist[0:3])
			
		elif templist[0] in ['.',',','*','A','a','T','t','C','c','G','g','N','n']:
			if len(templist)>1:                           ## 判断是不是最后一条reads的比对信息
				if templist[1] in ['+','-','$']:
					if templist[1]=="+":                  ## ,+2ta 插入
						ins_num = re.search("\d+",("").join(templist)).group()
						ins_seq = templist[2+len(ins_num):2+len(ins_num)+int(ins_num)]
						sign = ("").join(templist[0:2+len(ins_num)+int(ins_num)])
						del(templist[0:len(sign)])
						
					elif templist[1]=="-":                ## ,-2Aa 缺失
						del_num = re.search("\d+",("").join(templist)).group()
						del_seq = templist[2+len(del_num):2+len(del_num)+int(del_num)]
						# del_seq = re.search("[ACGTNacgtn]+",("").join(templist)).group()
						sign = ("").join(templist[0:2+len(del_num)+int(del_num)])
						# sign = templist[0]+re.search("\-\d+[ACGTNacgtn]+",("").join(templist)).group()
						del(templist[0:len(sign)])
						
					elif templist[1]=="$":                ## .$  末尾碱基
						sign = ("").join(templist[0:2])
						del(templist[0:2])
				else:                                     ## . , 正常比对到正负链或者 SNP
					sign=templist[0]
					del(templist[0])
			else:                                         ## 最后一条reads的碱基的比对情况，只能是 . , *
				sign=templist[0]
				del(templist[0])
		# elif templist[0] in ['A','a','T','t','C','c','G','g','N','n']: ## SNP
			# sign=templist[0]
			# del(templist[0])
		else:
			print "ERROR",chrom,ref_coordinate,mapinfo,index
		mapinfo_list[index]=sign
	return mapinfo_list
	
def pre_pileup(pileup,dictvcf):
	dictpileup={}
	dictSNVs={}
	dictallreads={}
	for line in open(pileup):
		chrom,ref_coordinate,refbase,depth,mapinfo,quality,reads_coordinate,reads=line.strip().split('\t')
		if len(dictvcf)>0:
			if ref_coordinate in dictvcf:
				reads_list=reads.split(",")
				reads_coordinate_list=reads_coordinate.split(",")
				mapinfo_list = process_pileup(chrom,ref_coordinate,refbase,depth,mapinfo,quality,reads_coordinate,reads)
				
				for n in range(int(depth)):
					readName=reads_list[n]
					dictallreads[readName]=""
					readsPos=reads_coordinate_list[n]
					mapInfo=mapinfo_list[n]
					if readName not in dictpileup:
						dictpileup[readName]={}
					if ref_coordinate not in dictpileup[readName]:
						dictpileup[readName][ref_coordinate]={}
					if mapInfo in [".",","]:
						mapInfo=refbase
					dictpileup[readName][ref_coordinate]=mapInfo.upper()
					
			# ALL pileup
			# if chrom not in dictpileup:
				# dictpileup[chrom]={}
			# if ref_coordinate not in dictpileup[chrom]:
				# dictpileup[chrom][ref_coordinate]={}
			# for n in range(int(depth)):
				# readName=reads_list[n]
				# readsPos=reads_coordinate_list[n]
				# mapInfo=mapinfo_list[n]
			# dictpileup[chrom][ref_coordinate][readName]=[readsPos,mapInfo]
		
		#SNVs
		'''
		if chrom not in dictSNVs:
			dictSNVs[chrom]={}
		if ref_coordinate not in dictSNVs[chrom]:
			dictSNVs[chrom][ref_coordinate]={}
		if re.search("[ACGTNacgtn]",mapinfo):
			for Match in mapinfo_list:
				if re.match("[ACGTNacgtn]",Match):##SNP
					if Match.upper() not in dictSNVs[chrom][ref_coordinate]:
						dictSNVs[chrom][ref_coordinate][Match.upper()]=0
					dictSNVs[chrom][ref_coordinate][Match.upper()]+=1
					
				indelInfo = re.search("[+-][0-9]+[ACGTNacgtn]+",Match):##indel
				if indelInfo:
					if indelInfo.group().upper() not in dictSNVs[chrom][ref_coordinate]:
						dictSNVs[chrom][ref_coordinate][indelInfo.group().upper()]=0
					dictSNVs[chrom][ref_coordinate][indelInfo.group().upper()]+=1
					
		'''
	return dictpileup,dictallreads

def getposinfo(dictvcf,dictpileup,pos_list,dictallreads,prefix,outdir): # modified by zhangsw
	'''
	fo=open("out.xls","w")
	fo.write("chrom	reads_name")
	for pos in pos_list:
		fo.write("\t"+str(pos))
	fo.write("\n")
	fo.write("chr6	type1")
	for pos in pos_list:
		fo.write("\t"+dictvcf[pos][0])
	fo.write("\n")
	fo.write("chr6	type2")
	for pos in pos_list:
		fo.write("\t"+dictvcf[pos][1])
	fo.write("\n")
	
	for readName in dictallreads:
		fo.write("chr6	"+readName)
		for pos in pos_list:
			if pos in dictpileup[readName]:
				fo.write("\t"+dictpileup[readName][pos])
			else:
				fo.write("\t"+"-")
		fo.write("\n")
	fo.close()
	'''
	## diff
	dictMapinfo={}
	for pos in pos_list:
		for ID in dictpileup:
			if pos in dictpileup[ID]:
				readbase=dictpileup[ID][pos]##dictpileup[readName][ref_coordinate]=mapInfo.upper()
				if ID not in dictMapinfo:
					dictMapinfo[ID]={}
					dictMapinfo[ID]["1"]=0
					dictMapinfo[ID]["2"]=0
				if readbase==dictvcf[pos][0]:
					dictMapinfo[ID]["1"]+=1
				if readbase==dictvcf[pos][1]:
					dictMapinfo[ID]["2"]+=1
	f1=open(outdir+"/"+prefix+".type1.list","w") # modified by zhangsw
	f2=open(outdir+"/"+prefix+".type2.list","w") # modified by zhangsw
	for name in dictMapinfo:
		if dictMapinfo[name]["1"]>dictMapinfo[name]["2"]:
			f1.write(name+"\n")
		if dictMapinfo[name]["1"]<dictMapinfo[name]["2"]:
			f2.write(name+"\n")
	f1.close()
	f2.close()
def main():
	opt = get_args()
	
	dictvcf,pos_list = read_vcf(opt.vcf)
	dictpileup,dictallreads = pre_pileup(opt.pileup,dictvcf)
	getposinfo(dictvcf,dictpileup,pos_list,dictallreads,opt.prefix,opt.outdir) # modified by zhangsw
	

if __name__ == '__main__':
	main()
