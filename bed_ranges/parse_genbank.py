#!/usr/bin/env python3
# coding:utf-8


from Bio import SeqIO
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="Parse GenBank file and convert it to bed format")

    parser.add_argument('input', type=str, help="GenBank format file to be parsed")
    parser.add_argument('output', type=str, help="Output prefix")
    parser.add_argument("--element", type=str, default=None,
                        help="Comma-separated list of element types to be parsed. Default is all.")
    parser.add_argument("--id", type=str, default=None,
                        help="reference name")
    parser.add_argument("--reverse", action="store_true",
                        help="Use reverse complement strand instead of forward strand")

    usage = '''
Usage:

python parse_genbank.3.py test.gb prefix

'''
    parser.epilog = usage
    return parser.parse_args()


def parse_region(region_str:str)->(int,int):
    # 解析输入的区间，格式为"start-end"
    start, end = None, None
    if region_str is not None:
        try:
            start, end = map(int, region_str.split('-'))
        except ValueError:
            raise ValueError("Invalid region format. Should be 'start-end'.")

    return start, end


def find_previous_feature(feature, features):
    index = features.index(feature)
    previous_feature = features[index-1]
    return previous_feature.type


def main():
    args = parse_args()

    # 读取GenBank文件
    records = SeqIO.parse(args.input, "genbank")

    # 解析所有指定类型的元素
    elements = set(args.element.split(',')) if args.element is not None else None

    GEN_START, GEN_END = None, None
    PREV_END = 0

    # 遍历每个记录
    out = open(args.output + ".bed", "w")
    outfa = open(args.output + ".fa", "w")
    for record in records:
        if args.id is not None:
            chrom = str(args.id)
        else:
            chrom = record.description.split(' ')[0]
            if "synthetic" in record.description.split(' '):
                chrom = str(record.name)
        # 输出fasta文件
        outfa.write(">" + chrom + "\n")
        outfa.write(str(record.seq) + "\n")
        # 这里用一个列表来保存记录的特征，方便后续的上下游区间处理
        features = list(record.features)
        # 按照位置排序，这样才能正确地找到上下游区间
        features.sort(key=lambda x: x.location.start.position)
        print(features)

        for feature in record.features:
            if feature.type == 'source':
                GEN_START = feature.location.start.position
                GEN_END = feature.location.end.position
                # print(GEN_START, GEN_END)
                continue

            # 判断是否为指定的元素类型
            if elements is not None and feature.type not in elements:
                continue

            # 解析特征的位置和注释信息
            start = feature.location.start.position
            end = feature.location.end.position
            strand = feature.location.strand
            name = feature.type
            if args.reverse:
                strand = -1 if strand == 1 else 1
            # print(start, end, strand, name)

            # 判断这个区间是不是和上个区间相连接
            if start - PREV_END > 0:
                prev_name = find_previous_feature(feature, features) + '-' + name
                bed_info = [chrom, PREV_END, start, prev_name, ".", strand]
                out.write('\t'.join(map(str, bed_info)) + '\n')
                PREV_END = start
                bed_info = [chrom, start, end, name, ".", strand]
                out.write('\t'.join(map(str, bed_info)) + '\n')
            elif start == PREV_END:
                # 将位置和注释信息输出为bed格式
                # bed_info = [record.id, start, end, name, ".", strand]
                bed_info = [chrom, start, end, name, ".", strand]
                out.write('\t'.join(map(str, bed_info)) + '\n')

            if features.index(feature) == len(features)-1 and end < GEN_END:
                name = name + '-end'
                bed_info = [chrom, end, GEN_END, name, ".", strand]
                out.write('\t'.join(map(str, bed_info)) + '\n')

            PREV_END = end

    out.close()
    outfa.close()
    return None


if __name__ == "__main__":
    main()
    #parse_region("1-100")
