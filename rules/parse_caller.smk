## Normalize each VCF file for truvari usage

global SAMPLES
global REF_DICT
global REFV
global MERGE_CALLERS
global DETECT_CALLERS

def find_ref(wildcards):
    return REF_DICT[REFV]

rule norm_caller_vcf:
    input:
        vcf="results/{sample}/{sample}.{caller}.vcf.gz",
    output:
        insdel="results/{sample}/{sample}.{caller}.insdel.vcf.gz",
    params:
        ref = find_ref
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod truvari/4.3.1
        if [ {wildcards.caller} == longcallD ]
        then
            bcftools norm --multiallelics - --output-type v {input.vcf}| bcftools view -i "(SVTYPE=='INS'||SVTYPE=='DEL')&FILTER=='PASS'" -O v - | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
            tabix -p vcf {output.insdel}
            
        elif [ {wildcards.caller} == pav ]
        then
            python {PIPELINE_DIR}/scripts/Pav2SV.py {input.vcf} {wildcards.sample} | bcftools norm --multiallelics - --output-type v /dev/stdin | python {PIPELINE_DIR}/scripts/resolve.py /dev/stdin pav {params.ref} |  bcftools norm --check-ref s --fasta-ref {params.ref} -N -m-any > results/{wildcards.sample}/pav.tmp.vcf
            bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v results/{wildcards.sample}/pav.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
            tabix -p vcf {output.insdel}

        elif [ {wildcards.caller} == dipcall ]
        then
            bcftools norm --multiallelics - --output-type v {input.vcf} | python /net/eichler/vol28/projects/medical_reference/nobackups/Scripts/MedRef/parsers/Dipcall2SV.py /dev/stdin {params.ref} | bcftools sort -o /dev/stdout -O z - > {output.insdel}
            tabix -p vcf {output.insdel}
            
        elif [ {wildcards.caller} == hapdiff ]
        then
            bcftools norm --multiallelics - --output-type v {input.vcf} | python {PIPELINE_DIR}/scripts/resolve.py /dev/stdin svimasm | bcftools norm --check-ref s --fasta-ref {params.ref} -N -m-any > results/{wildcards.sample}/hapdiff.tmp.vcf
            bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v results/{wildcards.sample}/hapdiff.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
            tabix -p vcf {output.insdel}

        else
            bcftools norm --multiallelics - --output-type v {input.vcf} | python {PIPELINE_DIR}/scripts/resolve.py /dev/stdin {wildcards.caller} {params.ref} |  bcftools norm --check-ref s --fasta-ref {params.ref} -N -m-any | bcftools annotate -x 'INFO/AF,INFO/STRAND' > results/{wildcards.sample}/{wildcards.caller}.tmp.vcf
            bcftools view -i "(SVTYPE=='INS'||SVTYPE=='DEL')&FILTER=='PASS'" -O v results/{wildcards.sample}/{wildcards.caller}.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}           
            tabix -p vcf {output.insdel}
            
        fi
        
        rm results/{wildcards.sample}/{wildcards.caller}.tmp.vcf
        """

rule parse_caller:
    input:
        insdel= expand("results/{sample}/{sample}.{caller}.insdel.vcf.gz", sample=SAMPLES.index, caller=MERGE_CALLERS.split(',')),
    # message:
    #     "Caller normalization complete"
    output:
        caller_list='results/{sample}/caller_vcf_list.txt',
        summary='results/{sample}/{sample}.caller_summary.txt'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        import gzip

        fout = open(output.caller_list,'w')
        for a_vcf in input.insdel:
            print(a_vcf,file=fout)
        fout.close()

        with open(output.summary,'w') as fout:
            print('caller\tINS\tDEL\tTOTAL',file=fout)
            for caller,vcf in zip(MERGE_CALLERS.split(','),input.insdel):
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