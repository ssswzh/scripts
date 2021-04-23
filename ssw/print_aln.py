#!/usr/bin/env python

import pysam
from pyfaidx import Fasta
from optparse import OptionParser

VERSION = 0.1

def print_aln(ref_aln, read_aln, aln, length, rs, qs):
  num, remainder = divmod(len(ref_aln), length)
  rp, qp = rs, qs
  for i in range(num):
    start, end = length * i, length * i + length
    tmp = ref_aln[start: end]
    print("Target: %8d    %s    %d" % (rp, tmp, rp + length - tmp.count('-') - 1))
    rp = rp + length - tmp.count('-')
    print("                    %s" % aln[start: end])
    tmp = read_aln[start: end]
    print("Query:  %8d    %s    %d" % (qp, tmp, qp + length - tmp.count('-') - 1))
    qp = qp + length - tmp.count('-')
    print()
  if remainder:
    start = length * num
    tmp = ref_aln[start: ]
    print("Target: %8d    %s    %d" % (rp, tmp, rp + len(tmp) - tmp.count('-')))
    print("                    %s" % aln[start:])
    tmp = read_aln[start:]
    print("Query:  %8d    %s    %d" % (qp, tmp, qp + len(tmp) - tmp.count('-')))
    print()

def process(options):
  ref = Fasta(options.ref, rebuild=False)
  bam = pysam.AlignmentFile(options.bam, "rb")
  for aln in bam:
    if aln.is_unmapped: continue
    chrm, ts, te = aln.reference_name, aln.reference_start, aln.reference_end
    qs, qe = aln.query_alignment_start, aln.query_alignment_end
    strand = '-' if aln.is_reverse else '+'
    cigartuples = aln.cigartuples
    print("target_name: %s\nquery_name: %s" % (chrm, aln.query_name))
    print("strand: %c\ttarget_begin: %d\ttarget_end: %d\tquery_begin: %d\tquery_end: %d\n" % (strand, ts, te, qs, qe))
    ref_seq = ref[chrm][ts: te + 1].seq
    read_seq = aln.query_alignment_sequence
    ref_tmp, read_tmp, aln_tmp = [], [], []
    ref_pos, read_pos = 0, 0
    for cigar in cigartuples:
      letter, length = cigar[:]
      if letter == 0:
        for i in range(length):
          if ref_seq[ref_pos + i] == read_seq[read_pos + i]:  aln_tmp.append("|")
          else: aln_tmp.append("*")
          ref_tmp.append(ref_seq[ref_pos + i])
          read_tmp.append(read_seq[read_pos + i])
        ref_pos, read_pos = ref_pos + length, read_pos + length
      elif letter == 1:
        aln_tmp.append(" " * length)
        ref_tmp.append("-" * length)
        read_tmp.append(read_seq[read_pos: read_pos + length])
        read_pos += length
      elif letter == 2:
        aln_tmp.append(" " * length)
        ref_tmp.append(ref_seq[ref_pos: ref_pos + length])
        read_tmp.append("-" * length)
        ref_pos += length
    print_aln("".join(ref_tmp), "".join(read_tmp), "".join(aln_tmp), int(options.length), ts, qs)

def parse_command():
  usage = "Convert sam/bam file to BLAST-like alignment\n\npython print_aln.py -r ref.fa -i input.bam"
  parser = OptionParser(usage=usage, version=VERSION)
  parser.add_option("-r", dest="ref", help="input reference file")
  parser.add_option("-i", dest="bam", help="input bam file")
  parser.add_option("-l", dest="length", help="line length in BLAST-like output [80]", default=80)
  return parser.parse_args()

def main():
  (options, args) = parse_command()
  if not options.ref: print("Please input reference file!"), exit(1)
  if not options.bam: print("Please input bam file!"), exit(1)
  process(options)

if __name__ == "__main__":
  main()

