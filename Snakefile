import os
import sys

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

include: "rules/caller_detect.smk"
# include: "rules/parse_caller.smk"
# include: "rules/intra_sample_collapse.smk"
# include: "rules/intra_sample_stats.smk"

rule detect:
    input:
        expand('results/{sample}/caller_detect.done', sample=SAMPLES.index)
    message:
        "Caller detection complete for all samples"

rule sample_detect:
    input:
        get_caller_detect
    output:
        flag = touch('results/{sample}/caller_detect.done')




# rule normalize:
#         norm_vcf=rules.parse_caller.input,
#
# rule intra:
#     input:
#         persample=rules.intra_sample_collapse.input,
#
# rule stats:
#     input:
#         stats=rules.intra_sample_stats.input

