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

## Parameters for cross caller integration
CROSS_CALLER_PARAMS = config.get('CROSS_CALLER_PARAMS')
CROSS_CALLER_SETV = 'caller_merge_{0}'.format(len(MERGE_CALLERS.split(',')))

SAMPLES = pd.read_csv(ALN_TABLE_FILENAME, sep='\t', index_col=['SAMPLE'])
REF = config.get('REF')




include: "rules/caller_detect.smk"
include: "rules/parse_caller.smk"
include: "rules/intra_sample_collapse.smk"
# include: "rules/intra_sample_stats.smk"

rule all:
    input:
        expand("results/{sample}/all_done",sample=SAMPLES.index)

rule detect:
    input:
        expand('results/{sample}/caller_detect.done', sample=SAMPLES.index)
    message:
        "Caller detection complete for all samples"

rule norm:
    input:
        expand('results/{sample}/{sample}.caller_summary.txt',sample=SAMPLES.index)
    message:
        "Caller normalization complete for all samples"

rule cross_caller:
    input:
        expand('results/{sample}/{caller_merge}/truvari_collapsed.insdel.pav-supp.vcf.gz', sample=SAMPLES.index, caller_merge=CROSS_CALLER_SETV),
    message:
        "Cross caller integration complete"

rule gather_outputs_per_sample:
    input:
        rules.detect.input,
        rules.norm.input,
        rules.cross_caller.input
    output:
        flag = touch("results/{sample}/all_done")



