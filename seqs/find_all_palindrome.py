#!/export/software/conda/miniconda3/bin/python

# Usage: python2 find_all_palindrome.py file seq_column_number

import sys

def check_palin(word):
    for i in xrange(len(word)/2):
        if word[i] != word[-1*(i+1)]:
            return False
    return True

def all_palindromes(string, pos): # modified by Zhangsiwen, 20191025

    left,right=0,len(string)
    j=right
    results=[]

    while left < right-1:
        temp = string[left:j] #Time complexity O(k)
        j-=1

        if check_palin(temp):
            position=int(pos)+left # modified by Zhangsiwen, 20191025
            results.append(str(position)+":"+temp) # modified by Zhangsiwen, 20191025

        if j<left+2:
            left+=1
            j=right

    return list(set(results))


data_file = open(sys.argv[1])
out_file = open(sys.argv[1]+".palindrome", "w")
for line in data_file:
    seq = line.strip().split("\t")[int(sys.argv[2])]
    pos = line.strip().split("\t")[1]
    out_file.write(line.strip()+"\t"+str(all_palindromes(seq, pos))+"\n")
data_file.close()
out_file.close()

