#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2019/10/25
# @Author  : zhangsiwen
# @contact : zhang.siwen@puruijizhun.com
# @LastModified: 2022/09/07
# @ChangeLog
#     20191025, first version
#     20220907, re-write script, change input file type to fasta, add --cutoff --mode


import argparse


def GetArgs(): #20220907
    parser = argparse.ArgumentParser(description='Find palindrome in sequences', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--fasta', help='fasta file', action='store', dest='fasta', required=True)
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--cutoff', help='length cutoff for palindrome, default %(default)s', default=3, type=int, action='store', dest='cutoff', required=False)
    optional.add_argument('--mode', help='output all palindromes or merge ones fully overlapped, default %(default)s', default='merge', choices=['all','merge'], action='store', dest='mode', required=False)
    usage = '''Usage:
    python %(prog)s --fasta fasta --cutoff length_cutoff --mode merge > out.bed
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


def check_palin(word):
    for i in range(0,round(len(word)/2)):
        if word[i] != word[-1*(i+1)]:
            return False
    return True


def merge_regions(results):
    new_result = []
    results.sort(key = lambda x: x[0])
    prev = [0, 0, '']
    for pos in results:
        start, end, seq = int(pos[0]), int(pos[1]), pos[2]
        if prev[0] <= start < prev[1] and end <= prev[1]:
            continue
        else:
            if prev != [0, 0, '']:
                new_result.append([str(prev[0]), str(prev[1]), str(prev[2])])
            prev = [start, end, seq]
    new_result.append([str(prev[0]), str(prev[1]), str(prev[2])])
    return new_result


def all_palindromes(string, cutoff, mode): # modified by Zhangsiwen, 20191025
    left,right=0,len(string)
    j=right
    results=[]
    while left < right-1:
        temp = string[left:j] #Time complexity O(k)
        j-=1
        if check_palin(temp):
            position = left # modified by Zhangsiwen, 20191025
            if len(temp) > cutoff:
                results.append([position,position+len(temp),temp]) # 20220907
        if j<left+2:
            left+=1
            j=right
    if mode == 'merge':
        results = merge_regions(results)
    return results


def main():
    ''' '''
    args = GetArgs()
    data_file = open(args.fasta)
    for line in data_file: # 20220907
        if line.startswith('>'):
            chrom = line.strip('>').split(' ')[0]
            continue
        seq = line.strip()
        results = all_palindromes(seq, args.cutoff, args.mode)
        for r in results:
            print(chrom+"\t"+'\t'.join(r))
    data_file.close()


if __name__ == '__main__':
    main()
