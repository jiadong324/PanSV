import os
import sys
import gzip

import pandas as pd


#
# Global constants
#

PIPELINE_DIR = os.path.dirname(os.path.realpath(workflow.snakefile))


#
# Parameters
#

configfile: config.get('config_file', 'config.yaml')

ALN_TABLE_FILENAME = config.get('MANIFEST')
HEADERS_DICT = config.get('HEADERS')
REF_DICT = config.get('REF')
REFV = config.get('REFV')
CPR_DICT = config.get('CPRMASK')
SD_DICT = config.get('SEGDUP')
TR_DICT = config.get('TRMASK')
# GENE_DICT = config.get('GENE')
MERGE_CALLERS = config.get('MERGE_CALLERS')
DETECT_CALLERS = config.get('DETECT_CALLERS')


SAMPLES = pd.read_csv(ALN_TABLE_FILENAME, sep='\t', index_col=['SAMPLE'])
REF = config.get('REF')


def get_caller_detect(wildcards):
    vcf_list = []
    sample = wildcards.sample

    for caller in DETECT_CALLERS.split(','):
        vcf_list.append(f'results/{sample}/{sample}.{caller}.vcf.gz')

    return vcf_list

def get_caller_norm(wildcards):
    vcf_list = []
    sample = wildcards.sample

    for caller in DETECT_CALLERS.split(','):
        vcf_list.append(f'results/{sample}/{sample}.{caller}.insdel.vcf.gz')

    return vcf_list

include: "rules/caller_detect.smk"
include: "rules/parse_caller.smk"
# include: "rules/intra_sample_collapse.smk"
# include: "rules/intra_sample_stats.smk"

rule detect:
    input:
        # expand('results/{sample}/caller_norm.done', sample=SAMPLES.index)
        expand('results/{sample}/caller_vcf_list.txt', sample=SAMPLES.index),
        expand('results/{sample}/{sample}.caller_summary.txt', sample=SAMPLES.index)

    message:
        "Caller detection and normalization complete for all samples"

rule run_caller:
    input:
        get_caller_detect
    output:
        flag = touch('results/{sample}/caller_detect.done')

rule norm_caller:
    input:
        insdel = get_caller_norm
    output:
        caller_list='results/{sample}/caller_vcf_list.txt',
        summary='results/{sample}/{sample}.caller_summary.txt'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        fout = open(output.caller_list,'w')
        for a_vcf in input.insdel:
            print(a_vcf,file=fout)
        fout.close()

        with open(output.summary,'w') as fout:
            print('caller\tINS\tDEL\tTOTAL',file=fout)
            for caller, vcf in zip(DETECT_CALLERS.split(','),input.insdel):
                ins = 0
                dele = 0
                with gzip.open(vcf,'rt') as fin:
                    for line in fin:
                        if line.startswith('#'):
                            continue
                        info = line.rstrip('\n').split('\t')[7]
                        info_fields = dict(kv.split('=',1) for kv in info.split(';') if '=' in kv)
                        svtype = info_fields.get('SVTYPE')
                        if svtype == 'INS':
                            ins += 1
                        elif svtype == 'DEL':
                            dele += 1
                total = ins + dele
                print(f'{caller}\t{ins}\t{dele}\t{total}',file=fout)



# rule intra:
#     input:
#         persample=rules.intra_sample_collapse.input,
#
# rule stats:
#     input:
#         stats=rules.intra_sample_stats.input

