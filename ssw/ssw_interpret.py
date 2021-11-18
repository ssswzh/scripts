#!/usr/bin/env python

import sys

if len(sys.argv) != 2:
    print("python " + sys.argv[0] + " ssw_test_result")
    print("Result will be in : ssw_test_result.flat and ssw_test_result.summary")
    exit
 
# primer aln info
stat = {}
target, query, score, strand = None, None, 0, None
fp = open(sys.argv[1], "r")
for line in fp:
    '''
    target_name: HPV16|NC_001526.4|E7gene
    query_name: HPV16_1
    optimal_alignment_score: 40	strand: +	target_begin: 120	target_end: 139	query_begin: 1	query_end: 20

    Target:      120    TCCAGCTGGACAAGCAGAAC    139
                        ||||||||||||||||||||
    Query:         1    TCCAGCTGGACAAGCAGAAC    20
    '''
    if line.startswith("target_name"): # ref
        target = line.strip().split(" ")[-1] # refname
    elif line.startswith("query_name"): # primer
        query = line.strip().split(" ")[-1] # p1/p2
        name = "_".join(query.split("_")[:-1]) # primer_id
        if name not in stat:
            stat[name] = {}
        if target not in stat[name]:
            stat[name][target] = []
        stat[name][target].append(query)
    elif line.startswith("optimal_alignment_score"):
        score = line.strip().split("\t")[0].split(" ")[-1] 
        stat[name][target].append(score)
        for ele in line.strip().split("\t"):
            if ele.startswith("query_length"):
                stat[name][target].append(ele.split(" ")[-1])
            if ele.startswith("strand"):
                stat[name][target].append(ele.split(" ")[-1])
            if ele.startswith("target_begin"):
                stat[name][target].append(ele.split(" ")[-1])
            if ele.startswith("target_end"):
                stat[name][target].append(ele.split(" ")[-1])
    # stat[primer_id] = {ref1: [p1, p1score, p1strand, p1start, p1end, p1length, p2, p2score, p2strand, p2start, p2end, p2length]}
    # stat[primer_id] = {ref1: [ 0,    1   ,     2   ,    3    ,  4  ,     5   , 6 ,    7   ,     8   ,    9   ,  10  ,    11   ]}
    
fp.close()
flat = open(sys.argv[1] + ".flat", "w")
flat.write("Query\tTarget\tP1\tP1_score\tP1_strand\tP1_start\tP1_end\tP1_length\tP2\tP2_score\tP2_strand\tP2_start\tP2_end\tP2_length\tTotal_score\n")
summary = open(sys.argv[1] + ".summary", "w")
summary.write("Query\tTarget\tP1_score\tP1_length\tP1_start\tP1_end\tP2_score\tP2_length\tP2_start\tP2_end\tTotal_score\n")
# store refname with highest score 
re = {}
for query,value in stat.items():
    re[query] = []
    # stat[primer_id] = {ref1: [p1, p1score, p1strand, p1start, p1end, p1length, p2, p2score, p2strand, p2start, p2end, p2length]}
    for ref,info in value.items():
        sums = str(int(info[1]) + int(info[7]))
        info[1] = str(info[1])
        info[7] = str(info[7])
        if info[2] == "+" and info[8] == "-":
            if re[query] == []:
                re[query] = [ref, info[1], info[5], info[3], info[4], info[7], info[11], info[9], info[10], sums]
            else:
                if int(sums) > int(re[query][-1]):
                    re[query] = [ref, info[1], info[5], info[3], info[4], info[7], info[11], info[9], info[10], sums]
            flat.write("\t".join([str(query),str(ref),"\t".join(info),str(sums)]) + "\n")
    #re[query][1] = str(re[query][1])
    summary.write("\t".join([str(query),"\t".join(re[query])]) + "\n")

flat.close()
summary.close()

