#!/usr/bin/env python3

import re
import argparse


def GetArgs():
    parser = argparse.ArgumentParser(description='Compare full length read map info to splited read map info, give correctly mapped read number.')
    parser.add_argument('--full', dest='full', help='full-length read map info', required=True) 
    parser.add_argument('--split', dest='split', help='splited read map info', required=True) 
    parser.add_argument('--distance', dest='distance', default=100, type=int, help='the distance between full-length read mapped position and splited read mapped position, default 100')
    parser.add_argument('--out', dest='output', help='output file name', required=True) 
    args = parser.parse_args()
    return args


def checkKEYinfo(read, d, info):
    # read: read_name
    # d: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    chrom, chrstart, chrend = info
    if read in d:
        if chrom in d[read]:
            d[read][chrom].append([chrstart, chrend])
        else:
            d[read][chrom] = [[chrstart, chrend]]
    else:
        d[read] = {chrom: [[chrstart, chrend]]}


def AlnPosFlatten(dic, distance):
    # dic: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    # same as d in checkKEYinfo(read, d, info)
    newdic = {}
    for read in dic:
        newdic[read] = {}
        for chrom in dic[read]:
            l = dic[read][chrom]
            # [[63074325,63077244], [63077244,63083955], [63083960,63085364], [77146724,77147929], [107126268,107128311], [107128311,107131313], [107131313,107134585]]
            l = sorted(l, key=lambda x:x[0])
            newl = []
            start, end = l[0]
            for r in l[1:]:
                if start - distance <= r[0] <= end + distance:
                    if r[1] >= end + distance:
                        end = r[1]
                else:
                    newl.append([start, end])
                    start, end = r
            newl.append([start, end]) # the last region stored in variable start and end
            newdic[read][chrom] = newl
    return newdic


def AlnInfo(infile, distance):
    aln = {}
    for line in open(infile,"r"):
        if line.startswith("ReadID"):
            continue
        elem = line.strip().split("\t")
        read_name, chrom, chrstart, chrend, strand = elem[0], elem[4], int(elem[5]), int(elem[6]), elem[9]
        readstart, readend, aln_length, read_length = int(elem[1]), int(elem[2]), int(elem[3]), int(elem[7])
        read_info = [chrom, chrstart, chrend]
        if read_name.endswith("ccs"): # full length read aln
            checkKEYinfo(read_name, aln, read_info)
        else: # splited read aln
            repos = re.search(r"/[0-9]+$", read_name).span()[0]
            read_id = read_name[0:repos]
            checkKEYinfo(read_id, aln, read_info)
    # flatten list of positions
    aln = AlnPosFlatten(aln, distance)
    return aln


def AlnCompare(full_aln, split_aln, distance, outfile):
    # full_aln: {read_name:{chrom:[[pos1, pos2], [...]], chrom:[[pos1, pos2], [...]]}, read_name2:{}}
    # newl = [[63074325, 63085364], [77146724, 77147929], [107126268, 107134585]]
    # a = [[63074325, 63085354], [77146724, 77147919], [107126268, 107134785]]
    # b = [[63074325, 63085354], [107126268, 107134785]]
    out = open(outfile, "w")
    out.write("ReadID\tFull-length Read Align Info\tSplited Read Align Info\tDecision\n")
    full_key_number = len(full_aln.keys())
    split_key_number = len(split_aln.keys())
    correct = 0 # cover true chrom region
    half_correct = 0 # cover sub-region of true chrom region
    partial_correct = 0 # cover one or more of true chrom regions, also mapped to other chrom or other region
    incorrect = 0 # unmapped or none of the chrom or region correct
    for read in full_aln:
        countTRUEchr = 0 # count correct chr number
        countFALSEchr = 0 # count incorrect chr number
        countHALFchr = 0 # count half_correct chr number
        countPARTIALchr = 0 # count partial_correct chr number
        if read in split_aln: # split has the same read
            full_chrom = list(full_aln[read].keys())
            split_chrom = list(split_aln[read].keys())
            for chrom in full_aln[read]:
                if chrom in split_aln[read]: # split has the same chrom of read
                    if len(full_aln[read][chrom]) == len(split_aln[read][chrom]): # full_aln in chrom has same region number compared to split_aln
                        countTRUEregion = 0 # count correct region number
                        countFALSEregion = 0 # count incorrect region number
                        countHALFregion = 0 # count half_correct region number
                        for i in range(0, len(full_aln[read][chrom])):
                            full_region = full_aln[read][chrom][i]
                            split_region = split_aln[read][chrom][i]
                            if abs(full_region[0]-split_region[0])<distance and abs(full_region[1]-split_region[1])<distance: # start and end difference less than distance
                                countTRUEregion = countTRUEregion + 1 # add one true chrom region
                            elif full_region[0]-distance < split_region[0] and split_region[1] < full_region[1]+distance: # start and end inside full_region, count as half_correct
                                countHALFregion = countHALFregion + 1 # add one half_correct chrom region
                            else:
                                countFALSEregion = countFALSEregion + 1 # add one incorrect chrom region
                        if countTRUEregion == len(full_aln[read][chrom]): # all region are correct
                            countTRUEchr = countTRUEchr + 1 
                        elif countFALSEregion == len(full_aln[read][chrom]): # all region are incorrect
                            countFALSEchr = countFALSEchr + 1 
                        elif countHALFregion != 0 and countFALSEregion == 0: # some region are half_correct, none is incorrect
                            countHALFchr = countHALFchr + 1
                        else: # no sub-region, only correct or incorrect regions
                            countPARTIALchr = countPARTIALchr + 1 
                    elif len(full_aln[read][chrom]) > len(split_aln[read][chrom]): # full_aln has more region than split_aln
                        countPARTIALchr = countPARTIALchr + 1 
                    elif len(full_aln[read][chrom]) < len(split_aln[read][chrom]): # split_aln has more region than full_aln
                        if full_aln[read][chrom][0][0]-distance < split_aln[read][chrom][0][0] and split_aln[read][chrom][-1][1] < full_aln[read][chrom][-1][1]+distance: # all split region inside full_aln region
                            countHALFchr = countHALFchr + 1
                        else: # some split region outside full_aln region
                            countFALSEchr = countFALSEchr + 1 
                else: # lack of a chrom
                    countFALSEchr = countFALSEchr + 1 # count incorrect chr number
            # final decision for read
            decision = ""
            if countTRUEchr != 0 and countFALSEchr == 0 and countHALFchr == 0 and countPARTIALchr == 0 and full_chrom == split_chrom: # full_aln and split_aln of read are same
                correct = correct + 1
                decision = "correct"
            elif countHALFchr != 0 and countFALSEchr == 0 and countPARTIALchr == 0 and full_chrom == split_chrom: # only half_correct, no incorrect or partial_correct, and aligned to same chroms
                half_correct = half_correct + 1
                decision = "half_correct"
            elif countTRUEchr == 0 and countHALFchr == 0 and countPARTIALchr == 0: # only incorrect
                incorrect = incorrect + 1
                decision = "incorrect"
            else: # all other conditions are partial_correct
                partial_correct = partial_correct + 1
                decision = "partial_correct"
            # print read info 
            out.write(read + "\t" + str(full_aln[read]) + "\t" + str(split_aln[read]) + "\t" + decision + "\n")
        else: # split unmapped
            incorrect = incorrect + 1
            out.write(read + "\t" + str(full_aln[read]) + "\t-\tincorrect\n")
    out.close()
    summary = open(outfile + ".summary", "w")
    summary.write("Batch\tTotalReads\tSplitedReads\tCorrect\tHalfCorrect\tPartialCorrect\tIncorrect\n")
    summary.write(outfile + "\t" + str(full_key_number) + "\t" + str(split_key_number) + "\t" + str(correct) + "\t" + str(half_correct) + "\t" + str(partial_correct) + "\t" + str(incorrect) + "\n")
    summary.write(outfile + "\t" + str(full_key_number/full_key_number*100) + "%\t" + str(split_key_number/full_key_number*100) + "%\t" + str(correct/full_key_number*100) + "%\t" + str(half_correct/full_key_number*100) + "%\t" + str(partial_correct/full_key_number*100) + "%\t" + str(incorrect/full_key_number*100) + "%\n")
    summary.close()
    

def main():
    args = GetArgs()
    full_aln = AlnInfo(args.full, args.distance)
    split_aln = AlnInfo(args.split, args.distance)
    AlnCompare(full_aln, split_aln, args.distance, args.output)


if __name__ == '__main__':
    main()


'''
countTRUEchr	countHALFchr	countPARTIALchr	countFALSEchr	list(full_aln[read].keys())==list(split_aln[read].keys())	decision
+	+	+	+	TRUE	partial_correct
+	+	+	-	TRUE	partial_correct
+	+	-	+	TRUE	partial_correct
+	+	-	-	TRUE	half_correct
+	-	+	+	TRUE	partial_correct
+	-	+	-	TRUE	partial_correct
+	-	-	+	TRUE	partial_correct
+	-	-	-	TRUE	correct
-	+	+	+	TRUE	partial_correct
-	+	+	-	TRUE	partial_correct
-	+	-	+	TRUE	partial_correct
-	+	-	-	TRUE	half_correct
-	-	+	+	TRUE	partial_correct
-	-	+	-	TRUE	partial_correct
-	-	-	+	TRUE	incorrect
-	-	-	-	TRUE	incorrect
+	+	+	+	FALSE	partial_correct
+	+	+	-	FALSE	partial_correct
+	+	-	+	FALSE	partial_correct
+	+	-	-	FALSE	partial_correct
+	-	+	+	FALSE	partial_correct
+	-	+	-	FALSE	partial_correct
+	-	-	+	FALSE	partial_correct
+	-	-	-	FALSE	partial_correct
-	+	+	+	FALSE	partial_correct
-	+	+	-	FALSE	partial_correct
-	+	-	+	FALSE	partial_correct
-	+	-	-	FALSE	partial_correct
-	-	+	+	FALSE	partial_correct
-	-	+	-	FALSE	partial_correct
-	-	-	+	FALSE	incorrect
-	-	-	-	FALSE	incorrect

'''

