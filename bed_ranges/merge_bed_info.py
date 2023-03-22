#!/usr/bin/env python3
# coding: utf-8
# author: zhang.siwen
# data: 2023.03.15


import sys
import os
import argparse
from collections import OrderedDict


def GetArgs():
    parser = argparse.ArgumentParser(description='Merge bed file info by coordinates (first three columns)', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--bed', help='input bed file', action='store', dest='bed', required=True)
    required.add_argument('--out', help='output file name with suffix and path', action='store', dest='out', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--field', help='column number for merging, give numbers separated by comma"," OR a range separated by colon ":",\nDEFAULT "%(default)s"\ne.g. 6 | 4,5,8,10 | 4:8', default='4:-1', action='store', dest='field', required=False)
    # optional.add_argument('--start', help='column number for starting merge, DEFAULT "%(default)s"', default=4, action='store', dest='start', type=int, required=False)
    # optional.add_argument('--end', help='column number for ending merge, DEFAULT "%(default)s"', default=-1, action='store', dest='end', type=int, required=False)
    args = parser.parse_args()
    return args


def ReadBed(bedfile, field='4:-1'):
    # merge field
    if ':' in field:
        if '-' in field:
            total_field = int(os.popen("tail -1 %s | awk -F'\t' '{print NF}'" % bedfile).readline().strip("\n"))
            split_field = [int(i) for i in field.strip().split(':')]
            field = []
            for i in split_field:
                if i > 0:
                    field.append(i)
                else:
                    field.append(total_field+i+1)
        else:
            field = [int(i) for i in field.strip().split(':')]
        field = list(range(field[0], field[1]+1))
    elif ',' in field:
        field = [int(i) for i in field.strip().split(',')]
    else:
        try:
            field = int(field)
        except:
            sys.exit('Please give a single number OR numbers separated by "," OR a range with two numbers separated by ":".')
    # header
    header = os.popen("head -1 %s " % bedfile).readline().strip("\n")
    # header
    if header.startswith('#'):
        header = header.strip().split('\t')
        header = header[0:3] + ['length'] + [header[i-1] for i in field]
    else:
        header = None
    # read bed
    bed = open(bedfile)
    bed_content = OrderedDict()
    for line in bed:
        if line.startswith('#'):
            continue
        # content
        element = line.strip().split('\t')
        mchr, mstart, mend = element[0:3]
        key = ':'.join([mchr, mstart, mend])
        if key in bed_content:
            for idx in field:
                bed_content[key][idx].append(element[idx-1])
        else:
            bed_content[key] = OrderedDict()
            bed_content[key]['length'] = str(int(mend) - int(mstart))
            for idx in field:
                bed_content[key][idx] = [element[idx-1]]
    bed.close()
    return bed_content, header


def SummarizeBed(bed_content, outfile='out', header=None):
    out = open(outfile, 'w')
    if header != None:
        out.write('\t'.join(header) + '\n')
    for key in bed_content.keys():
        length = bed_content[key]['length']
        del bed_content[key]['length']
        info_field = []
        for value in bed_content[key].values():
            content = set(value)
            if len(content) > 1:
                value = '|'.join(list(value))
            else:
                value = '|'.join(list(content))
            info_field.append(value)
        all_field = list(key.split(':')) + [length] + info_field
        out.write('\t'.join(all_field) + '\n')
    out.close()
    return


def Main():
    args = GetArgs()
    bed_content, header = ReadBed(args.bed, field=args.field)
    SummarizeBed(bed_content, outfile=args.out, header=header)


if __name__ == '__main__':
    Main()
