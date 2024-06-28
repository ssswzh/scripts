import argparse
import gzip
import concurrent.futures
import json


def get_args():
    parser = argparse.ArgumentParser(description="Calculate basic statistics from FASTQ files")
    parser.add_argument('-i', required=True, help="Input raw FASTQ file (single-end) or forward raw FASTQ file (paired-end)")
    parser.add_argument('-I', help="Input reverse raw FASTQ file (paired-end)")
    parser.add_argument('-c', required=True, help="Input clean FASTQ file (single-end) or forward clean FASTQ file (paired-end)")
    parser.add_argument('-C', help="Input reverse clean FASTQ file (paired-end)")
    parser.add_argument('-o', required=True, help="Output TSV file base name")
    parser.add_argument('--sample', required=True, help="Sample name")
    parser.add_argument('-j', help="Input JSON file with summary statistics")

    parser.epilog = '''
USAGE:
单端测序数据（无JSON文件）
    python3 %(prog)s -i raw_R1.fastq -c clean_R1.fastq -o output --sample RD-0R-0663-01-1
双端测序数据（无JSON文件）
    python3 %(prog)s -i raw_R1.fastq -I raw_R2.fastq -c clean_R1.fastq -C clean_R2.fastq -o output --sample RD-0R-0663-01-1
单端测序数据（有JSON文件）
    python3 %(prog)s -i raw_R1.fastq -c clean_R1.fastq -o output --sample RD-0R-0663-01-1 -j summary.json
双端测序数据（有JSON文件）
    python3 %(prog)s -i raw_R1.fastq -I raw_R2.fastq -c clean_R1.fastq -C clean_R2.fastq -o output --sample RD-0R-0663-01-1 -j summary.json

''' 
    return parser.parse_args()


def parse_fastq(file):
    if file.endswith('.gz'):
        open_func = lambda f: gzip.open(f, 'rt', encoding='utf-8')
    else:
        open_func = lambda f: open(f, 'rt', encoding='utf-8')

    with open_func(file) as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline().strip()
            qual = f.readline().strip()
            yield seq, qual


def calculate_n_stats(filenames, total_bases):
    total_n = 0

    def process_file(filename):
        local_n = 0
        for seq, _ in parse_fastq(filename):
            local_n += seq.count('N')
        return local_n

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(process_file, filenames)

    total_n = sum(results)
    n_ppm = round(total_n / total_bases * 1e6, 2) if total_bases > 0 else 0

    return n_ppm


def calculate_stats_from_files(filenames):
    total_reads = 0
    total_bases = 0
    total_gc = 0
    total_n = 0
    q20_bases = 0
    q30_bases = 0
    total_length = 0

    def process_file(filename):
        local_reads = 0
        local_bases = 0
        local_gc = 0
        local_n = 0
        local_q20 = 0
        local_q30 = 0
        local_length = 0

        for seq, qual in parse_fastq(filename):
            local_reads += 1
            seq_length = len(seq)
            local_bases += seq_length
            local_gc += seq.count('G') + seq.count('C')
            local_n += seq.count('N')
            local_q20 += sum(1 for q in qual if ord(q) - 33 >= 20)
            local_q30 += sum(1 for q in qual if ord(q) - 33 >= 30)
            local_length += seq_length

        return (local_reads, local_bases, local_gc, local_n, local_q20, local_q30, local_length)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = executor.map(process_file, filenames)

    for result in results:
        local_reads, local_bases, local_gc, local_n, local_q20, local_q30, local_length = result
        total_reads += local_reads
        total_bases += local_bases
        total_gc += local_gc
        total_n += local_n
        q20_bases += local_q20
        q30_bases += local_q30
        total_length += local_length

    avg_length = total_length // total_reads if total_reads > 0 else 0
    q20_percent = round(q20_bases / total_bases * 100, 2) if total_bases > 0 else 0
    q30_percent = round(q30_bases / total_bases * 100, 2) if total_bases > 0 else 0
    gc_percent = round(total_gc / total_bases * 100, 2) if total_bases > 0 else 0
    n_ppm = round(total_n / total_bases * 1e6, 2) if total_bases > 0 else 0

    return avg_length, total_reads, total_bases, q20_percent, q30_percent, gc_percent, n_ppm


def calculate_stats_from_json(json_file, before_filtering=True):
    with open(json_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    if before_filtering:
        summary = data['summary']['before_filtering']
    else:
        summary = data['summary']['after_filtering']
    
    avg_length = (summary['read1_mean_length'] + summary.get('read2_mean_length', 0)) // 2
    total_reads = summary['total_reads']
    total_bases = summary['total_bases']
    q20_percent = summary['q20_rate'] * 100
    q30_percent = summary['q30_rate'] * 100
    gc_percent = summary['gc_content'] * 100
    n_ppm = None  # N(ppm) will be calculated separately from FASTQ files
    
    return avg_length, total_reads, total_bases, q20_percent, q30_percent, gc_percent, n_ppm


def write_stats(filename, sample_name, stats):
    with open(filename, 'w', encoding='utf-8') as out_file:
        out_file.write("Sample\tLength(nt)\t#Reads\t#Bases\tQ20(%)\tQ30(%)\tGC(%)\tN(ppm)\n")
        out_file.write(f"{sample_name}\t" + "\t".join(map(str, stats)) + "\n")


def write_comparison(filename, sample_name, raw_stats, clean_stats):
    raw_reads, raw_bases = raw_stats[1], raw_stats[2]
    clean_reads, clean_bases = clean_stats[1], clean_stats[2]
    reads_ratio = (clean_reads / raw_reads * 100) if raw_reads > 0 else 0
    bases_ratio = (clean_bases / raw_bases * 100) if raw_bases > 0 else 0
    
    with open(filename, 'w', encoding='utf-8') as out_file:
        out_file.write("Sample\tRaw Reads\tClean Reads\tClean Reads Ratio (%)\tRaw Bases\tClean Bases\tClean Bases Ratio (%)\n")
        out_file.write(f"{sample_name}\t{raw_reads}\t{clean_reads}\t{reads_ratio:.2f}\t{raw_bases}\t{clean_bases}\t{bases_ratio:.2f}\n")


def main():
    args = get_args()
    
    raw_files = [args.i]
    if args.I:
        raw_files.append(args.I)
    
    clean_files = [args.c]
    if args.C:
        clean_files.append(args.C)
    
    if args.j:
        raw_stats = calculate_stats_from_json(args.j, before_filtering=True)
        clean_stats = calculate_stats_from_json(args.j, before_filtering=False)
        
        raw_n_ppm = calculate_n_stats(raw_files, raw_stats[2])
        clean_n_ppm = calculate_n_stats(clean_files, clean_stats[2])
        
        raw_stats = raw_stats[:-1] + (raw_n_ppm,)
        clean_stats = clean_stats[:-1] + (clean_n_ppm,)
    else:
        raw_stats = calculate_stats_from_files(raw_files)
        clean_stats = calculate_stats_from_files(clean_files)
    
    raw_output = args.o + ".raw.tsv"
    clean_output = args.o + ".clean.tsv"
    ratio_output = args.o + ".ratio.tsv"
    
    write_stats(raw_output, args.sample, raw_stats)
    write_stats(clean_output, args.sample, clean_stats)
    write_comparison(ratio_output, args.sample, raw_stats, clean_stats)


if __name__ == "__main__":
    main()
