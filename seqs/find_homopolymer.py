#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/09/07
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/09/07
# @Reference:
#     https://www.biostars.org/p/379454/#379465
# @ChangeLog
#     20220907, first version


import sys

if len(sys.argv) < 3:
    sys.exit('''Usage:
    python find_homopolymer.py fasta length_cutoff > out.bed
''')

inFile = open(sys.argv[1],'r')
l = int(sys.argv[2])
chrome = base = start = end = -1

for line in inFile:
    if line[0] == '>':
        if end - start >= l:
            print('\t'.join(map(str,[chrome,start,end,base*(end-start)])))
        chrome = line.strip().split()[0][1:]
        base = -1
        start = end = 0
    else:
        for b in line.strip():
            if b != base:
                if end - start >= l:
                    print('\t'.join(map(str,[chrome,start,end,base*(end-start)])))
                start = end
                base = b
            end += 1

if end - start >= l:
    print('\t'.join(map(str,[chrome,start,end,base*(end-start)])))

