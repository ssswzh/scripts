#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen@grandomics.com
# History:
#   Jan 7 2020, first version, to solve medaka_variant result combines two adjacent single position into one result, separate the result

import argparse

def GetArgs():
    des = "Compare two vcf files, generate TP, FP, FN, precision, recall, Fscore."
    parser = argparse.ArgumentParser(description = des)
    parser.add_argument('--base', dest='base', help='baseline vcf file', required=True) 
    parser.add_argument('--test', dest='test', help='test vcf file', required=True) 
    parser.add_argument('--prefix', dest='prefix', help='output prefix', required=True)
    parser.add_argument('--outdir', dest='outdir', help='output path', required=True) 
    args = parser.parse_args()
    return args

def VCFFormatConversion(vcffile, prefix):
    vcf = open(vcffile)
    out = open(prefix+".plain", "w")
    records = []
    for line in vcf:
        if not line.startswith("#"):
            ele = line.strip().split("\t")
            item = '_'.join([ele[i] for i in (0,1,3,4)])
            out.write(item+"\n")
            records.append(item)
    vcf.close()
    out.close()
    return records

def CompareVCFs(baseline, test, prefix, outdir):
    tp_out = open(outdir+"/"+prefix+".tp", "w")
    tp = set(baseline) & set(test)
    tp_out.write("\n".join(list(tp)))
    tp_out.close()
    fp_out = open(outdir+"/"+prefix+".fp", "w")
    fp = set(test) - set(baseline)
    fp_out.write("\n".join(list(fp)))
    fp_out.close()
    fn_out = open(outdir+"/"+prefix+".fn", "w")
    fn = set(baseline) - set(test)
    fn_out.write("\n".join(list(fn)))
    fn_out.close()
    precision = float(len(tp)/(len(tp)+len(fp)))
    recall = float(len(tp)/(len(tp)+len(fn)))
    Fscore = 2*precision*recall/(precision+recall)
    stats = open(outdir+"/"+prefix+".stats", "w")
    stats.write("\t".join(["Sample", "TP", "FP", "FN", "Precision", "Recall", "Fscore"])+"\n")
    stats.write("\t".join([prefix, str(len(tp)), str(len(fp)), str(len(fn)), str(precision), str(recall), str(Fscore)])+"\n")
    stats.close()

def main():
    args = GetArgs()
    basefile = args.base
    testfile = args.test
    prefix = args.prefix
    outdir = args.outdir
    base_format = VCFFormatConversion(basefile, basefile+".format")
    test_format = VCFFormatConversion(testfile, testfile+".format")
    CompareVCFs(base_format, test_format, prefix, outdir)

if __name__=="__main__":
    main()
