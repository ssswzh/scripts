#!/usr/bin/env python3
# coding: utf-8
# Author: zhangsiwen, May 31 2019
# Reference:
#   https://github.com/bmvdgeijn/WASP/blob/master/mapping/snptable.py
# History:
#   20201217, first version


import pysam
import argparse


# codes for CIGAR string
BAM_CMATCH     = 0   # M - match/mismatch to ref M
BAM_CINS       = 1   # I - insertion in read relative to ref
BAM_CDEL       = 2   # D - deletion in read relative to ref
BAM_CREF_SKIP  = 3   # N - skipped region from reference (e.g. intron)
BAM_CSOFT_CLIP = 4   # S - soft clipping (clipped sequence present in seq)
BAM_CHARD_CLIP = 5   # H - hard clipping (clipped sequence NOT present in seq)
BAM_CPAD       = 6   # P - padding (silent deletion from padded reference)
BAM_CEQUAL     = 7   # = - sequence match
BAM_CDIFF      = 8   # X - sequence mismatch



def GetArgs():
    parser = argparse.ArgumentParser(description='Parse CIGAR in bam, give mismatch, ins, del in each read.')
    parser.add_argument('--bam', dest='bam', help='bam file', required=True) 
    parser.add_argument('--out', dest='out', help='out file prefix', required=True)
    args = parser.parse_args()
    return args
    
    

def ParseCigar(bam):
    # read bam file
    samfile = pysam.Samfile(bam,"rb")

    # snp and indel positions
    positions = []
    # alt list
    alt = []
    # ref list
    ref = []
    # read id list
    readid = []
    
    # check each aln
    for record in samfile:
    
        # mapped reads
        if not record.is_unmapped:
        
            # init pos
            read_start = 0
            read_end = 0
            genome_start = 0
            genome_end = 0
            
            # cigar tuples
            for cigar in record.cigartuples:
                op = cigar[0] # CIGAR 'operation'
                op_len  = cigar[1] # length of operation
                
                # match or mismatch to reference
                if (op == BAM_CMATCH) or (op == BAM_CEQUAL) or (op == BAM_CDIFF):
                    
                    # check each base in 'xxM'
                    for i in range(0,op_len):
                        read_start = read_end
                        read_end = read_start + 1
                        genome_start = genome_end
                        genome_end = genome_start + 1 
                        a = record.query_sequence[read_start:read_end].upper()
                        r = record.get_reference_sequence()[genome_start:genome_end].upper()
                        
                        # mismatch inside this so-called Match pair
                        if a != r:
                            positions.append(record.reference_name + "\t" + str(record.pos+genome_start+1))
                            alt.append(str(a))
                            ref.append(str(r))
                            readid.append(record.query_name)
                
                # insert in read relative to reference
                elif op == BAM_CINS:
                    read_start = read_end
                    read_end = read_start + op_len
                    positions.append(record.reference_name + "\t" + str(record.pos+genome_start+1))
                    alt.append(str(record.query_sequence[read_start:read_end]))
                    ref.append("-")
                    readid.append(record.query_name)
                
                # deletion in read relative to reference
                elif op == BAM_CDEL:
                    genome_start = genome_end
                    genome_end   = genome_start + op_len
                    positions.append(record.reference_name + "\t" + str(record.pos+genome_start+1))
                    alt.append("-")
                    ref.append(str(record.get_reference_sequence()[genome_start:genome_end]))
                    readid.append(record.query_name)
                
                # section of skipped reference, such as intron
                elif op == BAM_CREF_SKIP:
                    genome_end = genome_end + op_len
                    genome_start = genome_end
                
                # this part of read skipped
                elif op == BAM_CSOFT_CLIP:
                    read_start = read_end
                    read_end = read_start + op_len
                
                # these bases not included in read or genome
                elif op == BAM_CHARD_CLIP:
                    pass
                
                # like an insert, likely only used in multiple-sequence
                # alignment where inserts may be of different lengths
                # in different seqs
                elif op == BAM_CPAD:
                    read_start += read_end
                    read_end = read_start + op_len

    return positions, ref, alt, readid



def WriteParsed(positions, ref, alt, readid, outfile):
    outParse = open(outfile, "w")
    for i in range(0,len(positions)):
        outParse.write("\t".join([positions[i], ref[i], alt[i], readid[i]]))
        outParse.write("\n")
    outParse.close()



def main():
    args = GetArgs()
    positions, ref, alt, readid = ParseCigar(args.bam)
    WriteParsed(positions, ref, alt, readid, args.out)



if __name__=="__main__":
    main()

