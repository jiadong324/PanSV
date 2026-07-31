

global SAMPLES
global DETECT_CALLERS
global REFV
global REF_DICT

def find_ref(wildcards):
    return REF_DICT[REFV]

def get_h1(wildcards):
    # manifest_df = pd.read_csv(MANIFEST[wildcards.ref],sep='\t',index_col=['NAME'])
    return SAMPLES.at[wildcards.sample, 'HAP1']

def get_h2(wildcards):
    # manifest_df = pd.read_csv(MANIFEST[wildcards.ref],sep='\t',index_col=['NAME'])
    return SAMPLES.at[wildcards.sample, 'HAP2']
def find_bam(wildcards):
    return SAMPLES.at[wildcards.sample, REFV]


rule sawfish_discover:
    input:
        bam = find_bam,
    output:
        bcf = 'results/{sample}/sawfish_disc/candidate.sv.bcf'
    params:
        ref = find_ref
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod sawfish/0.12.4
        sawfish discover --threads {threads} --ref {params.ref} --bam {input.bam} --output-dir $( dirname {output.bcf} ) --clobber
        """

rule sawfish_call:
    input:
        bcf = 'results/{sample}/sawfish_disc/candidate.sv.bcf'
    output:
        vcf = 'results/{sample}/{sample}.sawfish.vcf.gz'
    resources:
        mem = 50,
        hrs = 24,
        disk_free = 1,
    threads: 4
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod sawfish/0.12.4
        sawfish joint-call --threads {threads} --sample $( dirname {input.bcf} ) --clobber --output-dir $( dirname {output.vcf} )
        """


rule sniffles:
     input:
        bam = find_bam,
     output:
         vcf="results/{sample}/{sample}.sniffles.vcf.gz"
     params:
        ref = find_ref
     log:
         "log/results/{sample}.sniffles.log",
     resources:
         mem=10,
         hrs=24,
         disk_free=1,
     threads: 4
     shell:
         """
         source /etc/profile.d/modules.sh
         module load modules modules-init modules-gs/prod modules-eichler/prod sniffles/2.2
         sniffles -i {input.bam} --reference {params.ref} --output-rnames -v results/{wildcards.sample}/{wildcards.sample}.sniffles.vcf -t {threads}
         bcftools sort -Oz -o {output.vcf} results/{wildcards.sample}/{wildcards.sample}.sniffles.vcf
         tabix -p vcf {output.vcf}
         """

rule delly:
    input:
        bam=find_bam,
    output:
        vcf = "results/{sample}/{sample}.delly.bcf"
    params:
        ref=find_ref
    log:
        "log/results/{sample}.cutesv.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod delly/1.2.6
        delly lr -y pb -o {output.vcf} -g {params.ref} {input.bam}
        """

rule delly_bcf:
    input:
        bcf = rules.delly.output.vcf
    output:
        vcf = "results/{sample}/{sample}.delly.vcf.gz"
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.12.0
        bcftools view {input.bcf} | bgzip -c > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {input.bcf}
        """

rule cuteSV:
    input:
        bam = find_bam,
    output:
        vcf = "results/{sample}/{sample}.cutesv.vcf.gz"
    params:
        ref = find_ref
    log:
        "log/results/{sample}.cutesv.log",
    resources:
        mem=10,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod cuteSV/2.1.0
        cuteSV -t {threads} --write_old_sigs --genotype -l 50 --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 {input.bam} {params.ref} results/{wildcards.sample}/{wildcards.sample}.cutesv.vcf results/{wildcards.sample}
        bcftools sort -Oz -o {output.vcf} results/{wildcards.sample}/{wildcards.sample}.cutesv.vcf
        """


rule pbsv_discover:
    input:
        bam = find_bam,
    output:
        sig = "results/{sample}/{sample}.svsig.gz",
    params:
        ref = find_ref
    log:
        "log/results/{sample}.pbsv.log",
    resources:
        mem=50,
        hrs=24,
        disk_free=1,
    threads: 1
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod pbconda/202403
        pbsv discover --hifi --tandem-repeats {params.ref} {input.bam} {output.sig}
        """


rule pbsv_call:
    input:
        sig = "results/{sample}/{sample}.svsig.gz",
    output:
        vcf= "results/{sample}/{sample}.pbsv.vcf.gz"
    params:
        ref = find_ref
    log:
        "log/results/{sample}.pbsv.log",
    resources:
        mem=50,
        hrs=24
    threads: 4
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod pbconda/202403
        pbsv call --hifi -m 50 -j {threads} {params.ref} {input.sig} results/{wildcards.sample}/{wildcards.sample}.pbsv.vcf
        rm {input.sig}
        bcftools sort -Oz -o {output.vcf} results/{wildcards.sample}/{wildcards.sample}.pbsv.vcf
        tabix -p vcf {output.vcf}
        """

rule longcallD:
    input:
        bam=find_bam,
    output:
        vcf="results/{sample}/{sample}.longcallD.vcf.gz"
    params:
        ref=find_ref
    resources:
        mem=lambda wildcards, attempt: attempt * 32,
        hrs=24,
    threads: 6
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod longcallD/0.0.11
        longcallD call --hifi --min-sv-len 50 -n {wildcards.sample} -t {threads} -o results/{wildcards.sample}/{wildcards.sample}.longcallD.vcf {params.ref} {input.bam}       
        """

rule debreak:
    input:
        bam = find_bam,
    output:
        vcf="results/{sample}/{sample}.debreak.vcf.gz"
    params:
        ref = find_ref
    resources:
        mem=lambda wildcards, attempt: attempt * 32,
        hrs=24,
    threads: 6
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod sniffles/2.2
        debreak --bam {input.bam} --outpath $( dirname {output.vcf} ) --rescue_large_ins --poa --ref {params.ref}
        bcftools sort -Oz -o {output.vcf} results/{wildcards.sample}/debreak.vcf
        tabix -p vcf {output.vcf}
        """

rule svision:
    input:
        bam=find_bam,
    output:
        vcf= "results/{sample}/{sample}.svision.vcf"
    params:
        ref = find_ref,
    log:
        "log/results/{sample}.svision.log",
    resources:
        mem=20,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod SVision/1.4
        SVision -o $( dirname {output.vcf} ) -s 10 -b {input.bam} -t {threads} -n {wildcards.sample} -g {params.ref} -m /net/eichler/vol28/projects/medical_reference/nobackups/svision_model/svision-cnn-model.ckpt
        mv result/{wildcards.sample}/{wildcards.sample}.svision.*.vcf {output.vcf}
        """

rule svisionpro:
    input:
        bam = find_bam,
    output:
        vcf="results/{sample}/{sample}.svisionpro.vcf"
    params:
        ref = find_ref,
    log:
        "log/results/{sample}.svision.log",
    resources:
        mem=30,
        hrs=24,
        disk_free=1,
    threads: 8
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod svision-pro/2.3
        SVision-pro --out_path $( dirname {output.vcf} ) --min_supp 10 --preset hifi --target_path {input.bam} --process_num {threads} --sample_name {wildcards.sample} --genome_path {params.ref} --model_path /net/eichler/vol28/projects/medical_reference/nobackups/SVision-pro/src/pre_process/model_liteunet_256_8_16_32_32_32.pth
        mv result/{wildcards.sample}/{wildcards.sample}.svision_pro_*.vcf {output.vcf}
        """

rule hapdiff:
    input:
        h1 = get_h1,
        h2 = get_h2,
    output:
        vcf="results/{sample}/hapdiff_phased.vcf.gz",
    params:
        ref = find_ref
    resources:
        mem=20,
        hrs=24,
        disk_free=1,
    threads: 6
    shell:
        """
        source /etc/profile.d/modules.sh
        module load modules modules-init modules-gs/prod modules-eichler/prod miniconda/4.12.0 hapdiff/0.9
        hapdiff.py --reference {params.ref} --pat {input.h1} --mat {input.h2} --out-dir $( dirname {output.vcf} ) -t {threads} --sample {wildcards.sample} --sv-size 50
        """

rule caller_detect:
    input:
        expand("results/{sample}/{sample}.{var_caller}.vcf.gz", sample=SAMPLES.index, var_caller=DETECT_CALLERS.split(',')),
        expand("results/{sample}/hapdiff_phased.vcf.gz", sample=SAMPLES.index)
    message:
        "Caller detection complete"