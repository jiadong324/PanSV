import pandas as pd
import pysam
import numpy as np
import json
import os
from intervaltree import IntervalTree
import lib

global SAMPLES
global HEADERS_DICT
global REF_DICT
global CPR_DICT
global REFV
global MERGE_CALLERS
global CROSS_CALLER_PARAMS

def find_header(wildcards):
    return HEADERS_DICT[REFV]

def find_ref(wildcards):
    return REF_DICT[REFV]
def find_cprmask(wildcards):
    return CPR_DICT[REFV]


rule bcftool_all:
    input:
        vcf="results/{sample}/caller_vcf_list.txt"
    output:
        outvcf="results/{sample}/caller_merge/insdel.tmp.vcf.gz"
    log:
        "log/results/{sample}.bcftool.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 5
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.12.0
        bcftools merge --thread {threads} --merge none --force-samples -O z -o {output.outvcf} --file-list {input.vcf}
        tabix -p vcf {output.outvcf}
        """

rule truvari:
    input:
        bcfvcf=rules.bcftool_all.output.outvcf
    output:
        removed="results/{sample}/caller_merge/removed.vcf.gz",
        collapse="results/{sample}/caller_merge/truvari_collapsed.insdel.vcf.gz"
    log:
        "log/results/{sample}.truvari.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod truvari/5.2.0
        truvari collapse -i {input.bcfvcf} -c {output.removed} --sizemin 50 --sizemax 100000 --gt het -k first --intra {CROSS_CALLER_PARAMS} | bcftools sort --max-mem 8G -O z -o {output.collapse}
        tabix -p vcf {output.collapse}
        """


rule parse_truvari_collapse:
    input:
        vcf=rules.truvari.output.collapse,
        headers=find_header,
        excl=find_cprmask
    params:
        callers=MERGE_CALLERS,
    output:
        filt_vcf="results/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf",
        pav_supp_vcf="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf",
        pav_read_supp_vcf="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-read-supp.vcf",
        filt_bed="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.bed.gz",
        pav_read_supp_bed="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-read-supp.bed.gz",
        pav_bed="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-only.bed.gz",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        python {PIPELINE_DIR}/scripts/FilterSV.py -i {input.vcf} -o $( dirname {output.filt_vcf} ) -c {params.callers} -e {input.excl} -a {input.headers} -n {wildcards.sample}
        """

rule index_filt_vcf:
    input:
        filt_vcf="results/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf",
        pav_supp_vcf="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf"
    output:
        gzip="results/{sample}/caller_merge/truvari_collapsed.insdel.filt.vcf.gz",
        gzip_pav="results/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf.gz",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.12.0
        
        bcftools sort -o /dev/stdout -O v {input.filt_vcf} | bgzip -c > {output.gzip}
        tabix -p vcf {output.gzip}
        rm {input.filt_vcf}

        bcftools sort -o /dev/stdout -O v {input.pav_supp_vcf} | bgzip -c > {output.gzip_pav}
        tabix -p vcf {output.gzip_pav}
        rm {input.pav_supp_vcf}    
        """


# rule intra_sample_collapse:
#     input:
#         expand("results/{sample}/caller_merge/truvari_collapsed.insdel.pav-supp.vcf.gz", sample=SAMPLES.index)
#
#     message:
#         "Cross caller integration complete"