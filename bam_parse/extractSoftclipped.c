// Original script written by Devon Ryan
// Modified by zhangsiwen, zhangsiwen@grandomics.com
// Logs:
//      June 20 2019: change read ID header, add read name, mapping chr, soft-clip length, mapping coordinates, cigar value
//      June 21 2019: added argument to choose output file format either FASTQ or FASTA
//      July 11 2019: added argument to choose if print cigar string
// Ref: https://www.biostars.org/p/125874/
//      https://github.com/dpryan79/SE-MEI/blob/master/extractSoftclipped.c
//
// git clone https://github.com/dpryan79/SE-MEI.git
// cd SE-MEI
// git submodule update --init
// make

#include "htslib/sam.h"
#include "htslib/bgzf.h"
#include <stdlib.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

void process_alignment(bam1_t *b, bam_hdr_t *hdr, FILE *of, int length, int fastq, int ci) {// modified by zhangsw June21 2019
    int32_t pos = b->core.pos, apos = 0;
    int i, j, op, op_len;
    uint32_t *cigar = bam_get_cigar(b);
    char int2char[16] = {'N', 'A', 'C', 'N', 'G', 'N', 'N', 'N',
                         'T', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

    char *prefix; // modified by zhangsw July11 2019
	if(fastq){ // modified by zhangsw July11 2019
	    prefix = "@"; 
	} else {
	    prefix = ">" ;
	}
	
    for(i=0; i<b->core.n_cigar; i++) {
        op=bam_cigar_op(cigar[i]);
        op_len = bam_cigar_oplen(cigar[i]);
        if(op != 4 || op_len < length) {
            if(bam_cigar_type(op)&2) pos += op_len;
            if(bam_cigar_type(op)&1) apos += op_len;
        } else {
            if(ci){ // modified by zhangsw July11 2019
			    fprintf(of, "%s%s chr:%s len:%i pos:%"PRId32" cigar:", prefix, bam_get_qname(b), hdr->target_name[b->core.tid], op_len, pos+1); // modified by zhangsw July11 2019
				for(j=0; j<b->core.n_cigar; j++) {
					fprintf(of, "%i%c", bam_cigar_oplen(cigar[j]),BAM_CIGAR_STR[bam_cigar_op(cigar[j])]);
				}
            } else { // modified by zhangsw July11 2019
                fprintf(of, "%s%s chr:%s len:%i pos:%"PRId32, prefix, bam_get_qname(b), hdr->target_name[b->core.tid], op_len, pos+1); // modified by zhangsw July11 2019
			}
			fprintf(of, "\n"); // modified by zhangsw June21 2019
			for(j=0; j < op_len; j++) {
				fprintf(of, "%c", int2char[bam_seqi(bam_get_seq(b), j+apos)]);
			}
            if(fastq){// modified by zhangsw June21 2019
                fprintf(of, "\n+\n");
                for(j=0; j < op_len; j++) {
                    fprintf(of, "%c", bam_get_qual(b)[j+apos]+33);
                }
            }
            fprintf(of, "\n");
        }
    }
}

void usage(char *prog) {
    printf("Usage: %s [OPTIONS] file.bam > reads.fasta.gz\n", prog);
    printf("\n\
This program iterates over alignments in a SAM/BAM file and looks for reads that\n\
have soft-clipped regions. It then takes those soft-clipped regions and outputs\n\
them in gzipped fasta format to the console (so that one can redirect to a new\n\
file).\n\
\n\
Options\n\
\n\
-q INT   The minimum MAPQ of an alignment to be considered for inclusion\n\
         (default 0).\n\
\n\
-l INT   The minimum length of the soft-clipped segment (default 10 bases).\n\
\n\
-fq BOOL Output in 1(FASTQ) or 0(FASTA) format (default 0).\n\
\n\
-ci BOOL Output 1(cigar) or 0(no cigar) format (default 0).\n\
\n\
");
} // modified by zhangsw June21 2019

int main(int argc, char *argv[]) {
    bam_hdr_t *header = NULL;
    bam1_t *b = NULL;
    htsFile *fp = NULL;
    FILE *of = NULL;
    int i, op, op_len, is_softclipped, length = 10, minMAPQ = 0, fastq = 0, ci = 0; // modified by zhangsw July11 2019
    uint32_t *cigar;

    if(argc==1) {
        usage(argv[0]);
        return 0;
    }
    for(i=1; i<argc; i++) {
        if(strcmp(argv[i], "-h") == 0) {
            usage(argv[0]);
            return 0;
        } else if(strcmp(argv[i], "-l") == 0) {
            length = atoi(argv[++i]);
        } else if(strcmp(argv[i], "-q") == 0) {
            minMAPQ = atoi(argv[++i]);
        } else if(strcmp(argv[i], "-fq") == 0) { // modified by zhangsw June21 2019
            fastq = atoi(argv[++i]); // modified by zhangsw June21 2019
        } else if(strcmp(argv[i], "-ci") == 0) { // modified by zhangsw July11 2019
            ci = atoi(argv[++i]); // modified by zhangsw July11 2019
		} else if(fp == NULL) {
            fp = sam_open(argv[i], "rb");
            header = sam_hdr_read(fp);
        } else {
            usage(argv[0]);
            return -1;
        }
    }

    if(fp == NULL) {
        usage(argv[0]);
        return -1;
    }
    of = popen("gzip -c", "w");

    b = bam_init1();
    while(sam_read1(fp, header, b) >= 0) {
        if(b->core.n_cigar == 1) continue;
        if(b->core.qual < minMAPQ) continue;

        //Is it even soft-clipped?
        is_softclipped = 0;
        cigar = bam_get_cigar(b);
        for(i=0; i<b->core.n_cigar; i++) {
            op=bam_cigar_op(cigar[i]);
            op_len = bam_cigar_oplen(cigar[i]);
            if(op == 4 && op_len >= length) is_softclipped = 1;
        }
        if(is_softclipped) process_alignment(b, header, of, length, fastq, ci); // modified by zhangsw July11 2019
    }

    pclose(of);
    bam_destroy1(b);
    sam_close(fp);

    return 0;
}
