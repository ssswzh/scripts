#!/usr/bin/env python

import sys

target, query, score, strand = None, None, 0, None

stat = {}
fp = open(sys.argv[1], "r")
for line in fp:
  if line.startswith("target_name"): # ref
    target = line.strip().split(" ")[-1] # refname
  elif line.startswith("query_name"): # primer
    query = line.strip().split(" ")[-1] # p1/p2
    name = "_".join(query.split("_")[:-1]) # primer_id
    if name not in stat: stat[name] = {}
    if target not in stat[name]:  stat[name][target] = []
    stat[name][target].append(query) 
  elif line.startswith("optimal_alignment_score"):
    score = line.strip().split("\t")[0].split(" ")[-1] 
    stat[name][target].append(score)
    for ele in line.strip().split("\t"):
      if ele.startswith("strand"):
        stat[name][target].append(ele.split(" ")[-1])
      if ele.startswith("target_begin"):
        stat[name][target].append(ele.split(" ")[-1])
      if ele.startswith("target_end"):
        stat[name][target].append(ele.split(" ")[-1])
  # stat[primer_id]={ref1: [p1, p1score, p1strand, p1start, p1end, p2, p2score, p2strand, p2start, p2end]}
  # stat[primer_id]={ref1: [ 0,    1   ,     2   ,    3   ,   4  , 5 ,    6   ,     7   ,    8   ,   9  ]}
      
fp.close()
re = {}
plain = open(sys.argv[2],"w")
summary = open(sys.argv[3],"w")
for key,value in stat.items():
  re[key] = []
  for k, v in value.items():
    sums = str(int(v[1]) + int(v[6]))
    v[1] = str(v[1])
    v[6] = str(v[6])
    if v[2] == "+" and v[7] == "-":
      if re[key] == []:
        re[key] = [k, sums, v[1], v[6], v[3], v[4], v[8], v[9]]
      else:
        if int(sums) > int(re[key][1]):
          re[key] = [k, sums, v[1], v[6], v[3], v[4], v[8], v[9]]
      plain.write("\t".join([str(key),str(k),"\t".join(v),str(sums)]) + "\n")
  #re[key][1] = str(re[key][1])
  summary.write("\t".join([str(key),"\t".join(re[key])]) + "\n")

plain.close()
summary.close()

