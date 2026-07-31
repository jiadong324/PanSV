## Normalize each VCF file for truvari usage

global SAMPLES
global REF_DICT
global REFV
global MERGE_CALLERS

def find_ref(wildcards):
    return REF_DICT[REFV]

rule norm_caller_vcf:
    input:
        vcf="results/{sample}/{sample}.{caller}.vcf",
    output:
        insdel="results/{sample}/{sample}.{caller}.insdel.vcf.gz",
    params:
        ref = find_ref
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.2.1",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python scripts/resolve.py /dev/stdin {wildcards.caller} {input.ref} |  bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any | bcftools annotate -x 'INFO/AF,INFO/STRAND' > results/{wildcards.sample}/{wildcards.caller}.tmp.vcf
        bcftools view -i "(SVTYPE=='INS'||SVTYPE=='DEL')&FILTER=='PASS'" -O v results/{wildcards.sample}/{wildcards.caller}.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}           
        tabix -p vcf {output.insdel}
        rm results/{wildcards.sample}/{wildcards.caller}.tmp.vcf
        """

rule parse_hapdiff:
    input:
        vcf = "results/{sample}/hapdiff_phased.vcf.gz",
    output:
        insdel="results/{sample}/{sample}.hapdiff.insdel.vcf.gz",
        # inv="{ref}/{sample}/{sample}.hapdiff.inv.vcf.gz",
    params:
        ref=find_ref
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "truvari/4.2.1",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python /net/eichler/vol28/projects/medical_reference/nobackups/Scripts/MedRef/parsers/resolve.py /dev/stdin svimasm | bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any > results/{wildcards.sample}/hapdiff.tmp.vcf
        bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v results/{wildcards.sample}/hapdiff.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        tabix -p vcf {output.insdel}
        rm results/{wildcards.sample}/hapdiff.tmp.vcf
        """

rule parse_pav:
    input:
        vcf='results/{sample}/pav_{sample}.vcf.gz',
    output:
        insdel="results/{sample}/{sample}.pav.insdel.vcf.gz",
        inv="results/{sample}/{sample}.pav.inv.vcf.gz",
    params:
        ref=find_ref
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "miniconda/4.12.0",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        python ../scripts/Pav2SV.py {input.vcf} {wildcards.sample} | bcftools norm --multiallelics - --output-type v /dev/stdin | python ../scripts/resolve.py /dev/stdin pav {input.ref} |  bcftools norm --check-ref s --fasta-ref {input.ref} -N -m-any > results/{wildcards.sample}/pav.tmp.vcf
        bcftools view -i "SVTYPE=='INS'||SVTYPE=='DEL'" -O v results/{wildcards.sample}/pav.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        bcftools view -i "SVTYPE=='INV'" -O v results/{wildcards.sample}/pav.tmp.vcf | bcftools sort -o /dev/stdout -O v - | bgzip -c > {output.insdel}
        tabix -p vcf {output.insdel}
        tabix -p vcf {output.inv}
        rm results/{wildcards.sample}/pav.tmp.vcf
        """

rule parse_dipcall:
    input:
        vcf='results/{sample}/{sample}.dip.vcf.gz',
    output:
        vcf='results/{sample}/{sample}.dip.insdel.vcf.gz'
    params:
        ref=find_ref
    envmodules:
        "modules",
        "modules-init",
        "modules-gs/prod",
        "modules-eichler/prod",
        "dipcall/0.3",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        bcftools norm --multiallelics - --output-type v {input.vcf} | python /net/eichler/vol28/projects/medical_reference/nobackups/Scripts/MedRef/parsers/Dipcall2SV.py /dev/stdin {input.ref} | bcftools sort -o /dev/stdout -O z - > {output.vcf}
        tabix -p vcf {output.vcf}
        """


rule caller_vcf_list:
    input:
        caller_vcf=expand("results/{sample}/{sample}.{caller}.insdel.vcf.gz", sample=SAMPLES.index, caller=MERGE_CALLERS.split(','))
    output:
        vcf_list = 'results/{sample}/vcf_list.txt'
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    run:
        fout = open(output.vcf_list, 'w')
        for line in input.caller_vcf:
            print(line, file=fout)


rule parse_caller:
    input:
        expand("results/{sample}/vcf_list.txt", sample=SAMPLES.index)
    message:
        "Caller normalization complete"