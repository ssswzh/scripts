#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import argparse

# Define the depth thresholds
DEPTH_THRESHOLDS = [1, 2, 5, 10, 50, 100, 150, 200]

# output table header
HEADER = [
    "Sample",
    "Total Reads",
    "Mapped Reads",
    "Mapped Reads Ratio (%)",
    "Unique Reads",
    "Unique Reads Ratio (%)",
    "Genome Bases",
    "Covered Bases",
    "Coverage (%)",
    "Mean Depth",
    "1X",
    "1~2X",
    "2~5X",
    "5~10X",
    "10~50X",
    "50~100X",
    "100~150X",
    "150~200X",
    "≥200X"
]


def get_args():
    # Get arguments
    parser = argparse.ArgumentParser(description="Summary bam info from picard result.")
    parser.add_argument("--alignsummary", "-a", help="Align summary file from CollectAlignmentSummaryMetrics.", required=True)
    parser.add_argument("--coverage", "-c", help="Coverage file from CollectTargetedPcrMetrics.", required=True)
    # parser.add_argument("--insert", "-c", help="InsertSize file from CollectInsertSizeMetrics.", required=True)
    parser.add_argument("--depth", "-d", help="depth file from samtools depth", required=True)
    parser.add_argument("--output", "-o", help="Out file, TSV format.", required=True)
    parser.add_argument("--sample", "-p", help="sample name", required=True)
    args = parser.parse_args()
    return args


def parse_metrics(metric_file, file_type):
    f = open(metric_file, "r")
    line = f.readline()
    metrics_dict = {}
    if file_type == "CollectAlignmentSummaryMetrics":
        while line:
            if line.startswith("CATEGORY"):
                header_list = line.strip().split("\t")
                f.readline()
                f.readline()
                values_list = f.readline().strip().split("\t")
                metrics_dict = dict(zip(header_list, values_list))
            line = f.readline()
    elif file_type == "CollectTargetedPcrMetrics" or file_type == "CollectInsertSizeMetrics":
        while line:
            if line.startswith("CUSTOM_AMPLICON_SET") or line.startswith("MEDIAN_INSERT_SIZE"):
                header_list = line.strip().split("\t")
                values_list = f.readline().strip().split("\t")
                metrics_dict = dict(zip(header_list, values_list))
            line = f.readline()
    return metrics_dict


def calculate_depth_percentages(file_path, thresholds):
    # Initialize counters for each depth range
    depth_counts = [0] * (len(thresholds) + 1)
    total_bases = 0
    with open(file_path, 'r') as file:
        for line in file:
            total_bases += 1
            depth = int(line.strip().split('\t')[2])

            if depth == thresholds[0]:
                depth_counts[0] += 1
            else:
                for i in range(1, len(thresholds)):
                    if thresholds[i-1] <= depth < thresholds[i]:
                        depth_counts[i] += 1
                        break
                else:
                    depth_counts[-1] += 1
    # Calculate percentages
    percentages = [round(count / total_bases * 100, 2) for count in depth_counts]
    return percentages


def main():
    args = get_args()
    output = args.output
    sample = args.sample

    alignsummary = parse_metrics(args.alignsummary, "CollectAlignmentSummaryMetrics")
    coverage = parse_metrics(args.coverage, "CollectTargetedPcrMetrics")
    # insert = parse_metrics(args.insert, "CollectInsertSizeMetrics")
    coverage_depth = calculate_depth_percentages(args.depth, DEPTH_THRESHOLDS)
    result = [
        str(args.sample),
        alignsummary["TOTAL_READS"],
        alignsummary["PF_READS_ALIGNED"],
        alignsummary["PCT_PF_READS_ALIGNED"],
        coverage["PF_UQ_READS_ALIGNED"],
        coverage["PCT_PF_UQ_READS_ALIGNED"],
        coverage["GENOME_SIZE"],
        coverage_depth[0],
        round(coverage_depth[0] / int(coverage["GENOME_SIZE"]) * 100, 2),
        round(float(coverage["MEAN_TARGET_COVERAGE"]))
        ] + coverage_depth
    out = open(args.output, 'w')
    for header, value in zip(HEADER, result):
        out.write(f"{header}\t{value}\n")
    out.close()


if __name__ == "__main__":
    main()
