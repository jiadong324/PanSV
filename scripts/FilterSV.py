#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 10/30/24
'''
import gzip
import argparse
import pandas as pd
import numpy as np
import json
import os
import pysam
from intervaltree import IntervalTree

def to_binary(n, bits):
    return ''.join(str(1 & int(n) >> i) for i in range(bits))

def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict

def get_overlaps(a_start, a_end, b_start, b_end):
    return max(0, min(a_end, b_end) - max(a_start, b_start))

def create_size_df(ins_length_info, ins_values, del_length_info, del_values, ref, sample):
    df_data_ins = pd.DataFrame()
    df_data_del = pd.DataFrame()

    num = 0
    for x, y in zip(ins_values,[ins_length_info[i] for i in ins_values]):
        num += 1
        df_data_ins.loc[num, "x"] = x
        df_data_ins.loc[num, "y"] = np.log10(y) if y != 0 else 0
        df_data_ins.loc[num, "type"] = "ins"
        df_data_ins.loc[num, "ref"] = ref
        df_data_ins.loc[num, "sample"] = sample
    for x, y in zip(del_values,[del_length_info[i] for i in del_values]):
        num += 1
        df_data_del.loc[num, "x"] = -x
        df_data_del.loc[num, "y"] = np.log10(y) if y != 0 else 0
        df_data_del.loc[num, "type"] = "del"
        df_data_del.loc[num, "ref"] = ref
        df_data_del.loc[num, "sample"] = sample
    df_data = pd.concat([df_data_ins, df_data_del])
    return df_data

def fetch_depth(depth_tabix, chrom, start, end):
    reads = []
    for line in depth_tabix.fetch(chrom, start, end):
        entries = line.strip().split('\t')
        reads.append(int(entries[3]))

    if not reads:
        return 0

    return int(np.mean(reads))


if __name__ == '__main__':


    AUTOSOMESXY = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", 'chrY']

    parser = argparse.ArgumentParser(
        description="Create intra-sample SV and add support evidence for each SV",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('-i', '--input', dest='vcf', help="Path to Truvari collapsed VCFs from multiple callers", required=True)
    parser.add_argument('-o', '--outdir', dest='outdir', help="Output directory for filtered results", required=True)
    parser.add_argument('-c', '--caller', dest='caller', help="Callers used for Truvari collapse", required=True)
    parser.add_argument('-e', '--excl', dest='excl', help="Path to excluded regions", required=False)
    parser.add_argument('-a', '--header', dest='header', help="Path to VCF header for filtered VCF file", required=True)
    parser.add_argument('-n', '--sample', dest='sample', help="Sample name", required=True)
    # parser.add_argument('-d', '--depth', dest='depth', help="Read depth of this sample", required=True)
    parser.add_argument('-q', '--qual', dest='qual', help="Path to the quality results predicted by BoostSV", required=False)


    read_based_callers = ['pbsv', 'cutesv', 'nanovar', 'sniffles', 'cutesv', 'delly', 'sawfish']
    other_callers = ['hapdiff', 'dipcall', 'pbsv', 'cutesv', 'nanovar', 'sniffles', 'cutesv', 'delly', 'sawfish']

    args = parser.parse_args()

    # depth_bed = pysam.Tabixfile(args.depth)

    qual_tbl = None
    if args.qual:
        qual_df = pd.read_csv(args.qual, sep='\t', header=[0])
        qual_tbl = qual_df.sort_values(by='QUAL', ascending=False)
        qual_tbl.drop_duplicates(subset='ID', keep="first", inplace=True)

        qual_tbl.set_index('ID', inplace=True)

    sample = args.sample
    caller_list = args.caller.split(',')

    outdir = args.outdir
    header_file = args.header

    ## Read excluded regions to intervaltree
    excl_dict = {}
    if args.excl:
        for line in open(args.excl):
            entries = line.strip().split('\t')
            chrom, start, end = entries[0], entries[1], entries[2]
            if chrom not in excl_dict:
                excl_dict[chrom] = IntervalTree()
            excl_dict[chrom][int(start): int(end)] = (int(start), int(end))

    insdel_dict = {'Total': {'INS': 0, 'DEL': 0}}
    caller_dict = {}
    sv_records = []
    svs_bed, pav_supp_bed, pav_read_supp_bed, pav_only_svs, qual_filt_bed = [], [], [], [], []
    filt_svs_vcf = []
    pav_supp_vcf, pav_read_supp_vcf = {}, []
    qual_filt_vcf = []
    alt_fa_out = open(f'{outdir}/truvari_collapsed.insdel.pav-supp.fa', 'w')
    alt_asm_unique_fa = open(f'{outdir}/truvari_collapsed.insdel.pav-only.fa', 'w')


    for line in gzip.open(args.vcf, 'rt'):
        if line.startswith('#'):
            continue
        entries = line.strip().split('\t')
        bin_supp = int(entries[9].split(':')[-1])
        supp_vec = to_binary(bin_supp, len(caller_list))
        supp = sum([int(ele) for ele in supp_vec])

        if entries[0] not in AUTOSOMESXY:
            continue
        if supp == 0:
            continue

        info_dict = parse_vcf_info_column(entries[7])
        svlen = abs(int(info_dict['SVLEN']))
        if svlen < 50 or svlen > 100000:
            continue

        svtype = info_dict['SVTYPE']
        end = int(entries[1]) + 1 if svtype == 'INS' else int(entries[1]) + svlen


        insdel_dict['Total'][svtype] += 1
        if f'SUPP={supp}' not in insdel_dict:
            insdel_dict[f'SUPP={supp}'] = {'INS': 0, 'DEL': 0}
        insdel_dict[f'SUPP={supp}'][svtype] += 1

        info_dict['SVLEN'] = -svlen if svtype == 'DEL' else svlen
        info_dict['END'] = end

        supp_caller = ''
        is_read_sv = False
        is_pav_sv = False
        is_supp_sv = False
        for i, val in enumerate(supp_vec):
            if val == '1':
                supp_caller += f'{caller_list[i]},'
                if caller_list[i] == 'pav':
                    is_pav_sv = True

                if caller_list[i] in other_callers:
                    is_supp_sv = True

                if caller_list[i] in read_based_callers:
                    is_read_sv = True

        if supp_caller[:-1] not in caller_dict:
            caller_dict[supp_caller[:-1]] = {'INS': 0, 'DEL': 0}
        caller_dict[supp_caller[:-1]][svtype] += 1


        info_dict['SUPP_CALLER'] = supp_caller[:-1]
        new_sv_id = f'{entries[0]}-{int(entries[1]) + 1}-{svtype}-{abs(svlen)}'

        info_dict['VARID'] = f'{new_sv_id}-{sample}-{bin_supp}'
        info_labels = ['VARID', 'SVLEN', 'END', 'SVTYPE', 'SUPP_CALLER']

        if is_pav_sv:
            info_labels.append("QRY_REGION" if "QRY_REGION" in info_dict else "TIG_REGION")

        info_list = []

        for key in info_labels:
            if key in ['QRY_REGION', 'TIG_REGION']:
                # print(info_dict)
                info_list.append('TIG_REGION={1}'.format(key, info_dict[key]))
                continue
            info_list.append('{0}={1}'.format(key, info_dict[key]))

        new_info_str = ';'.join(info_list)
        in_excl = '.'
        if args.excl and excl_dict[entries[0]].overlaps(int(entries[1]), end):
            in_excl = 'Exclude_Region'
            svs_bed.append([entries[0], entries[1], end, svlen, svtype, entries[2], new_sv_id, entries[9].split(':')[0],
                            supp_caller[:-1], in_excl])
            continue
        svs_bed.append([entries[0], entries[1], end, svlen, svtype, entries[2], new_sv_id, entries[9].split(':')[0], supp_caller[:-1], in_excl])


        ## PAV-only SVs
        if supp_caller[:-1] == 'pav':

            # depth = fetch_depth(depth_bed, entries[0], int(entries[1]) - 50, end + 50)
            phased_gt = entries[9].split(':')[0]
            sv_record = f'{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{bin_supp}\tPASS\t{new_info_str}\tGT:AL:APOS\t{phased_gt}:{abs(svlen)}:{entries[1]}'
            filt_svs_vcf.append(sv_record)
            pav_only_svs.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, entries[9].split(':')[0], supp_caller[:-1], sample])

            if svtype == 'INS':
                print(f'>{entries[2]}', file=alt_asm_unique_fa)
                print(f'{entries[4]}', file=alt_asm_unique_fa)
            if svtype == 'DEL':
                print(f'>{entries[2]}', file=alt_asm_unique_fa)
                print(f'{entries[3]}', file=alt_asm_unique_fa)

            if args.qual:
                sv_score = -1 if np.isnan(qual_tbl.at[new_sv_id, 'QUAL']) else qual_tbl.at[new_sv_id, 'QUAL']
                sig_reads = -1 if np.isnan(qual_tbl.at[new_sv_id, 'SIG_READS_NUM']) else qual_tbl.at[new_sv_id, 'SIG_READS_NUM']
                if sv_score >= 0.5:
                    qual_filt_bed.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, entries[9].split(':')[0], supp_caller[:-1], sample])
                    qual_info_str = f'{new_info_str};ML_QUAL={sv_score};SIG_READS_NUM={sig_reads}'
                    sv_record = f'{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{bin_supp}\tPASS\t{qual_info_str}\tGT:AL:APOS\t{phased_gt}:{abs(svlen)}:{entries[1]}'
                    qual_filt_vcf.append(sv_record)

        if is_pav_sv and is_supp_sv:
            # depth = fetch_depth(depth_bed, entries[0], int(entries[1]) - 50, end + 50)
            phased_gt = entries[9].split(':')[0]
            sv_record = f'{entries[0]}\t{entries[1]}\t{new_sv_id}\t{entries[3]}\t{entries[4]}\t{bin_supp}\tPASS\t{new_info_str}\tGT:AL:APOS\t{phased_gt}:{abs(svlen)}:{entries[1]}'

            if is_read_sv:
                pav_read_supp_vcf.append(sv_record)
                pav_read_supp_bed.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, phased_gt, supp_caller[:-1], sample])

            filt_svs_vcf.append(sv_record)

            pav_supp_bed.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, phased_gt, supp_caller[:-1], sample])

            if entries[1] not in pav_supp_vcf:
                pav_supp_vcf[entries[1]] = [sv_record, bin_supp]
            else:
                dup_rec_qual = pav_supp_vcf[entries[1]][1]
                if bin_supp > dup_rec_qual:
                    pav_supp_vcf[entries[1]] = [sv_record, bin_supp]

            if svtype == 'INS':
                print(f'>{new_sv_id}', file=alt_fa_out)
                print(f'{entries[4]}', file=alt_fa_out)
            if svtype == 'DEL':
                print(f'>{new_sv_id}', file=alt_fa_out)
                print(f'{entries[3]}', file=alt_fa_out)

            if args.qual:
                sv_score = -1 if np.isnan(qual_tbl.at[new_sv_id, 'QUAL']) else qual_tbl.at[new_sv_id, 'QUAL']
                sig_reads = -1 if np.isnan(qual_tbl.at[new_sv_id, 'SIG_READS_NUM']) else qual_tbl.at[new_sv_id, 'SIG_READS_NUM']
                if sv_score >= 0.5:
                    qual_filt_bed.append([entries[0], entries[1], end, svlen, svtype, new_sv_id, entries[9].split(':')[0],supp_caller[:-1], sample])
                    qual_info_str = f'{new_info_str};ML_QUAL={sv_score};SIG_READS_NUM={sig_reads}'
                    sv_record = f'{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{bin_supp}\tPASS\t{qual_info_str}\tGT:AL:APOS\t{phased_gt}:{abs(svlen)}:{entries[1]}'
                    qual_filt_vcf.append(sv_record)

    alt_asm_unique_fa.close()
    alt_fa_out.close()

    filt_sv_writer = open(f'{outdir}/truvari_collapsed.insdel.filt.vcf', 'w')
    pav_supp_vcf_writer = open(f'{outdir}/truvari_collapsed.insdel.pav-supp.vcf', 'w')
    pav_read_supp_vcf_writer = open(f'{outdir}/truvari_collapsed.insdel.pav-read-supp.vcf', 'w')

    for line in open(header_file):
        print(line.strip(), file=filt_sv_writer)
        print(line.strip(), file=pav_supp_vcf_writer)
        print(line.strip(), file=pav_read_supp_vcf_writer)

    if args.qual:
        print('##INFO=<ID=ML_QUAL,Number=1,Type=Float,Description="SV quality predicted by BoostSV">', file=filt_sv_writer)
        print('##INFO=<ID=SIG_READS_NUM,Number=1,Type=Integer,Description="Number of signature reads identified by BoostSV">', file=filt_sv_writer)

    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}', file=filt_sv_writer)
    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}', file=pav_supp_vcf_writer)
    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}', file=pav_read_supp_vcf_writer)

    for sv_line in filt_svs_vcf:
        print(sv_line, file=filt_sv_writer)

    filt_sv_writer.close()

    for sv_id, (sv_line, qual) in pav_supp_vcf.items():
        print(sv_line, file=pav_supp_vcf_writer)

    pav_supp_vcf_writer.close()

    for sv_line in pav_read_supp_vcf:
        print(sv_line, file=pav_read_supp_vcf_writer)

    pav_read_supp_vcf_writer.close()

    ## Save all collapsed SVs to BED file
    sorted_bed = sorted(svs_bed, key=lambda x: (x[0], int(x[1])))
    pd.DataFrame(sorted_bed, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'VCF_ID', 'GT', 'SUPP', 'EXCL']).to_csv(
        f'{outdir}/truvari_collapsed.insdel.bed.gz', sep='\t', header=False, index=False, compression='gzip')

    ## Save all collapsed SVs to BED file exclude masked regions
    filt_svs_bed = sorted(pav_supp_bed + pav_only_svs, key=lambda x: (x[0], int(x[1])))
    pd.DataFrame(filt_svs_bed, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'MERGE_SAMPLES']).to_csv(
        f'{outdir}/truvari_collapsed.insdel.filt.bed.gz', sep='\t', header=True, index=False, compression='gzip')

    ## Save PAV SVs supported by callers to BED file
    sorted_supp_bed = sorted(pav_supp_bed, key=lambda x: (x[0], int(x[1])))
    pd.DataFrame(sorted_supp_bed, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'MERGE_SAMPLES']).to_csv(
        f'{outdir}/truvari_collapsed.insdel.pav-supp.bed.gz', sep='\t', header=True, index=False, compression='gzip')

    ## Save PAV SVs supported by read callers to BED file
    sorted_read_supp_bed = sorted(pav_read_supp_bed, key=lambda x: (x[0], int(x[1])))
    pd.DataFrame(sorted_read_supp_bed,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'MERGE_SAMPLES']).to_csv(
        f'{outdir}/truvari_collapsed.insdel.pav-read-supp.bed.gz', sep='\t', header=True, index=False, compression='gzip')

    ## Save PAV only SVs to BED file
    sorted_pav_bed = sorted(pav_only_svs, key=lambda x: (x[0], int(x[1])))
    pd.DataFrame(pav_only_svs, columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'MERGE_SAMPLES']).to_csv(
        f'{outdir}/truvari_collapsed.insdel.pav-only.bed.gz', sep='\t', header=True, index=False, compression='gzip')

    ## Save SV counts by number of supporting callers
    with open(f'{outdir}/truvari_collapsed.insdel.stats.json', 'w') as f:
        json.dump(insdel_dict, f)

    ## Save SVs detected by different combination of callers
    with open(f'{outdir}/truvari_collapsed.insdel.callers.json', 'w') as f:
        json.dump(caller_dict, f)


    if args.qual:
        sorted_qual_bed = sorted(qual_filt_bed, key=lambda x: (x[0], int(x[1])))
        pd.DataFrame(sorted_qual_bed,columns=['#CHROM', 'POS', 'END', 'SVLEN', 'SVTYPE', 'ID', 'GT', 'SUPP', 'MERGE_SAMPLES']).to_csv(f'{outdir}/truvari_collapsed.insdel.boostsv_pass.bed.gz', sep='\t', header=True, index=False,
            compression='gzip')

        qual_filt_vcf_fout = open(f'{outdir}/truvari_collapsed.insdel.boostsv_pass.vcf', 'w')

        for line in open(header_file):
            print(line.strip(), file=qual_filt_vcf_fout)
        print('##INFO=<ID=ML_QUAL,Number=1,Type=Float,Description="SV quality predicted by BoostSV">', file=qual_filt_vcf_fout)
        print('##INFO=<ID=SIG_READS_NUM,Number=1,Type=Integer,Description="Number of signature reads identified by BoostSV">', file=qual_filt_vcf_fout)
        print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}', file=qual_filt_vcf_fout)

        for sv_line in qual_filt_vcf:
            print(sv_line, file=qual_filt_vcf_fout)

        qual_filt_vcf_fout.close()

    # size100_df = create_size_df(ins_length_info1, ins_values1, del_length_info1, del_values1, args.ref, sample)
    # size1k_df = create_size_df(ins_length_info2, ins_values2, del_length_info2, del_values2, args.ref, sample)
    #
    # size100_df.to_csv(f'{outdir}/truvari_collapsed.insdel.size100_1000.csv')
    # size1k_df.to_csv(f'{outdir}/truvari_collapsed.insdel.size1000_10000.csv')
