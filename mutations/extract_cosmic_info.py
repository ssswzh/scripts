#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2022/08/26
# @Author  : zhangsiwen
# @contact : zhang.siwen
# @LastModified: 2022/08/26
# @ChangeLog
#     20220826, first version


import argparse
from collections import OrderedDict


def GetArgs():
    parser = argparse.ArgumentParser(description='Extract COSMIC info to TAB-delimited file', formatter_class=argparse.RawTextHelpFormatter)
    parser._action_groups.pop()
    # required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument('--vcf', help='cosmic database vcf file, \noutput all records if --id and --ids not specified', action='store', dest='vcf', required=True)
    # optional arguments
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('--out', help='out file name', default=None, action='store', dest='out', required=False)
    optional.add_argument('--ids', help='give cosmic ids (both COSVxx and COSMxx are OK), \nif need to look for multiple ids, use comma "," to separate them, \ne.g. "COSV58737130,COSV58737076"', default=None, action='store', dest='ids', required=False)
    optional.add_argument('--idfile', help='a file containing all cosmic ids, one id per line', default=None, action='store', dest='idfile', required=False)
    optional.add_argument('--dedup', help='remove duplication records, keep only one record', default=False, action='store_true', dest='dedup', required=False)
    usage = '''Usage:
    %(prog)s --vcf CosmicCodingMuts.normal.vcf 
    %(prog)s --vcf CosmicCodingMuts.normal.vcf --out out.tsv --id COSV58737130,COSV58737076,COSM911918
    %(prog)s --vcf CosmicCodingMuts.normal.vcf --out out.tsv --idfile cosmic_id.txt
'''
    parser.epilog = usage
    args = parser.parse_args()
    return args


def ExtractCosmicDB(vcffile, cosmic_id=None, dedup=True):
    '''
    @func: read vcf file, store all information
    '''
    vcf = open(vcffile)
    result = OrderedDict()
    header_field = []
    for line in vcf:
        if line.startswith('#'):
            continue
        records = line.strip().split('\t')
        CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFOField = records[0:8]
        if header_field == []:
            header_field = [ele.split('=')[0] for ele in INFOField.split(';')]
        info = {ele.split('=')[0]:ele.split('=')[1] for ele in INFOField.split(';')}
        ids = [ID, info['LEGACY_ID']] 
        if cosmic_id:
            for i in ids:
                if i in cosmic_id:
                    if dedup:
                        if '_ENST' not in info['GENE']:
                            result[i] = [CHROM, POS, ID, REF, ALT] + list(info.values())
                    else:
                        result[i] = [CHROM, POS, ID, REF, ALT] + list(info.values())
        else:
            if dedup:
                if '_ENST' not in info['GENE']:
                    result[ID] = [CHROM, POS, ID, REF, ALT] + list(info.values())
            else:
                result[ID] = [CHROM, POS, ID, REF, ALT] + list(info.values())
    vcf.close()
    return result, header_field


def GetCosmicIDs(ids, idfile):
    if ids is not None:
        ids = ids.strip().split(',')
    elif idfile is not None:
        ids = []
        with open(idfile) as inputfile:
            for line in inputfile:
                ids.append(line.strip())
    return ids


def WriteOutput(result, header_field, outfile):
    if outfile:
        out = open(outfile, 'w')
        out.write('\t'.join(header_field) + '\n')
        for key, value in result.items():
            out.write('\t'.join(value) + '\n')
        out.close()
    else:
        for key, value in result.items():
            print('\t'.join(header_field))
            print('\t'.join(value))
    return
    

def main():
    ''' '''
    args = GetArgs()
    if args.ids or args.idfile:
        cosmic_id = GetCosmicIDs(args.ids, args.idfile)
        result, header_field = ExtractCosmicDB(args.vcf, cosmic_id, args.dedup)
    else:
        result, header_field = ExtractCosmicDB(args.vcf, args.ids, args.dedup)
    WriteOutput(result, header_field, args.out)


if __name__ == '__main__':
    main()
