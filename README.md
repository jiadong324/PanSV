# PanSV

The workflow used for detecting SVs from long-read genomes in population-scale and clinical.
It runs multiple different assembly- and read-based callers with a per sample level SV quality assessment.

Check [1KG_LongRead_SV](https://github.com/jiadong324/1KG_LongRead_SV) for the release of SVs detected by this pipeline from 1KG long-read genomes. 

## SV discovery

### Setup

1. Link file ```rundist``` to your working directory.
2. Modify the ```./config/config.yaml``` example file as needed and copy to your working directory.
3. Create the a manifest file of alignment file. The header should be ```NAME\tCHM13\tGRCh38```. Each caller will look for BAM files based on the ```REFV``` in ```config.yaml``` and it has to match the manifest header.

### Alignment
HiFi reads are aligned with [pbmm2](https://github.com/PacificBiosciences/pbmm2) v1.13.1 ‘--preset HiFi’. 
1KG-ONT reads are aligned with minimap2 in the NAPU pipeline used in previous publication by Jonas Gustafson et.al.
We used IB-ONT alignment directly from the publication for SV discovery and phasing.

Minimap2 v2.28 is used to align assembly to both references. The alignment pipeline is [here](https://github.com/mrvollger/asm-to-reference-alignment).

**NOTE:** This pipeline currently dose not support alignment.

### SV callers

To detect SVs for each genome, simply run ```./rundist detect 50``` or a dry-run with ```./rundist detect 50 -np```

This will create ```{sample}.{caller}.insdel.vcf``` and a summary table of SVs detected by different caller for each sample ```{sample}.caller_summary.txt```.

The ```{sample}.{caller}.insdel.vcf``` will be used to create a high-quality per genome SV callset.

The current pipeline supports the following callers. 

| Tool        | Input type | Version | Website                                       |
|-------------|------------|---------|-----------------------------------------------|
| PAV         | Assembly   | v2.3.4  | https://github.com/EichlerLab/pav             |
| Dipcall     | Assembly   | v0.3    | https://github.com/lh3/dipcall                |
| SVIM-ASM    | Assembly   | v0.9    | https://github.com/eldariont/svim-asm         |
| pbsv        | HiFi       | v2.9.0  | https://github.com/PacificBiosciences/pbsv    |
| sawfish     | HiFi       | v0.12.4 | https://github.com/PacificBiosciences/sawfish |
| sniffles    | HiFi, ONT  | v2.2    | https://github.com/fritzsedlazeck/sniffles    |
| delly       | HiFi, ONT  | v1.2.6  | https://github.com/dellytools/delly           |
| cutesv      | HiFi, ONT  | v2.1.0  | https://github.com/tjiangHIT/cuteSV           |
| Nanovar     | ONT        | v1.8.0  | https://github.com/benoukraflab/NanoVar       |
| debreak     | HiFi, ONT  | v1.2.0  | https://github.com/Maggi-Chen/DeBreak         |
| SVision     | HiFi, ONT  | v1.4    | https://github.com/xjtu-omics/SVision         |
| SVision-pro | HiFi, ONT  | v2.3    | https://github.com/songbowang125/SVision-pro  |
| LongcallD   | HiFi       | v0.0.11 | https://github.com/yangao07/longcallD         |

### TR genotyping

Tandem repeat catalogs for GRCh38 and T2T-CHM13 can be found [here](https://zenodo.org/records/13178746).


| Tool   | Input type | Version | Dataset             | Purpose        |
|--------|------------|---------|---------------------|----------------|
| TRGT   | HiFi       | v1.4.1  | HPRC, HGSVC, UW-ONT | Tandem repeats |
| vamos  | Assembly   | v2.1.5  | HPRC, HGSVC, UW-ONT | Tandem repeats |


### Per genome SV

For each genome, we prioritize the PAV calling results and identify SVs supported by at least one another caller with [Truvari (v5.2.0)](https://github.com/acenglish/truvari). 
The pipeline is ```rules/intra_sample_collapse.smk```. The output of this pipeline is a multi-caller integrated VCF used the PAV reported breakpoint, sv length, phased genotype etc.
We then run [BoostSV](https://github.com/jiadong324/BoostSV) on the multi-caller integrated VCF for each genome.

We also used the same annotation as [Logsdon et al. Nature 2025](https://www.nature.com/articles/s41586-025-09140-6) to exclude SVs inside complex regions, gaps, etc. 
Briefly, these regions include UCSC gaps and centromere on GRCh38. 
For T2T-CHM13, complex regions include centromere, acrocentric p-arms, satellite regions except for monomeric satellite.


**NOTE:** TR repeat genotypes are not included in the current persample SV calling output. More benchmarks have to be done for this to be added. 


## Cohort-level reference panel

**NOTE:** For the current SV reference panel created from HPRC/HGSVC phased assemblies, we did not do additional phasing of SVs and SNPs because they are called from same haplotype input.

### Callable regions
We first defined the callable regions for each genome. HGSVC/HPRC and UW-ONT callable regions were created by PAV based on the assembly to reference alignment. 
For IB-ONT genomes, we splitted the PMDV phased BAM into two read sets. The callable regions for each haplotype were created by merging each read set (>= 3 reads) into non–overlapping intervals with BEDtools merge ‘-d 500’. 

### Create callset

The non-redundant set integrated SVs from HPRC, HGSVC, UW-ONT and BI-ONT genomes with coverage >= 15x and read N50 >=15 kbp. 
Truvari (v5.2.0) was used to create the non-redundant SV set for both references with ‘--pctseq 0.90 –pctsize 0.90 –refdist 500 –keep common’. 
Briefly, different alleles at the same SV site were collapsed if they share minimum 90% sequence similarity and 90% allele size similarity. 
The calls with the highest quality predicted by BoostSV were used to represent each collapsed SV site. 
Moreover, this integration only considered INS/DEL ranging from 50bp to 100,000bp. 

We then filled the missing genotypes ‘./.’ with reference genotypes ‘0|0’, ‘0|.’ and ‘.|0’ for each genome based on callable regions. 
Note that the missing genotype in the final VCF only suggests there is no confident read or assembly alignments. 
For each integrated SV, we also kept the allele breakpoint position (FORMAT/APOS) and length (FORMAT/AL) from each sample. 
BCFtools (v1.16) plugin function fill-tags is used to calculate the statistics for each SV site, including allele frequency, minor allele frequency, etc (https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html). 