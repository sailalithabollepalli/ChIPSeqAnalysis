# Get sample names
with open('samples.txt') as f:
    SAMPLES= f.read().splitlines()

# Reference fasta file (hg19)
REF = "hg19.fa" 

# blacklisted peaks
BL= "wgEncodeDacMapabilityConsensusExcludable.bed"

rule all:
    input:
        expand("filtered_sorted_BL_reads/{sample}_sorted.dedup.q30.unique.sorted.BL.bam", sample=SAMPLES),
        expand("filtered_sorted_BL_reads/{sample}_sorted.dedup.q30.unique.sorted.BL.bam.readcounts", sample=SAMPLES),
        "peaks_macs/mergedSamples_macs_peaks.bed"

rule fastqc:
    input:
        "{sample}.fastq",
    output:
        "FastQC/{sample}_fastqc.zip",
        "FastQC/{sample}_fastqc.html"
    params: ""
    log:
        "logs/fastqc/{sample}.log"
    shell:
        "fastqc {input} -o ./FastQC"

rule multiqc:
    input:
        expand("FastQC/{sample}_fastqc.html",sample=SAMPLES),
        expand("FastQC/{sample}_fastqc.zip",sample=SAMPLES)
    output:
        "FastQC/multiqc_report.html"
    shell:
        "multiqc {input} -o ./FastQC"

rule bwa_map:
    input:
        ref=REF,
        file="{sample}.fastq"
    output:
        bam="mapped_reads/{sample}.bam"
    log:
        "logs/bwa_map/{sample}.log"
    shell:
        "bwa mem {input.ref} {input.file} | samtools view -Sb - > {output.bam}"

rule sort_bam:
    input:
        bamfile="mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}_sorted.bam"
    log:    
        "logs/samtools_sort/{sample}.log"
    params: 
        jobname = "{sample}"  
    shell:
        "samtools sort -T {output}.tmp -o {output} {input.bamfile} 2> {log}"

rule index_bam:
    input:
        "sorted_reads/{sample}_sorted.bam"
    output:
        "sorted_reads/{sample}_sorted.bam.bai"
    log:    
        "logs/samtools_index/{sample}.log"
    params: 
        jobname = "{sample}"  
    shell:
        "samtools index {input} 2> {log}"

rule samtools_dedup:
    input:
        "sorted_reads/{sample}_sorted.bam"
    output:
        "dedup_reads/{sample}_sorted.dedup.bam"
    shell:
        "samtools rmdup -s {input} {output}"

rule reads_filtering:
    input:
        "dedup_reads/{sample}_sorted.dedup.bam"
    output:
        "filtered_reads/{sample}_sorted.dedup.q30.unique.bam"
    shell:
        "samtools view -h -q 30 -F 4 -F 256 {input} | grep -v XA:Z | grep -v SA:Z | samtools view -b - > {output} "
                                        
rule samtools_sort2:
    input:
        bamfiltered="filtered_reads/{sample}_sorted.dedup.q30.unique.bam"
    output:
        "filtered_sorted_reads/{sample}_sorted.dedup.q30.unique.sorted.bam"
    log:    
        "logs/filtered_sortedreads/{sample}.log"
    params: 
        jobname = "{sample}"  
    shell:
        "samtools sort -T {output}.tmp -o {output} {input.bamfiltered} 2> {log}"

rule samtools_index2:
    input:
        "filtered_sorted_reads/{sample}_sorted.dedup.q30.unique.sorted.bam"
    output:
        "filtered_sorted_reads/{sample}_sorted.dedup.q30.unique.sorted.bam.bai"
    log:    
        "logs/filtered_sortedreads_index/{sample}.log"
    params: 
        jobname = "{sample}"  
    shell:
        "samtools index {input} 2> {log}"

rule blacklist_filtering:
    input:
        bam="filtered_sorted_reads/{sample}_sorted.dedup.q30.unique.sorted.bam",
        bed=BL
    output:
        "filtered_sorted_BL_reads/{sample}_sorted.dedup.q30.unique.sorted.BL.bam"
    params:
        jobname = "{sample}"  
    shell: 
        "bedtools intersect -v -a {input.bam} -b {input.bed} > {output}"

rule readcount:
    input:
        "filtered_sorted_BL_reads/{sample}_sorted.dedup.q30.unique.sorted.BL.bam"
    output:
        "filtered_sorted_BL_reads/{sample}_sorted.dedup.q30.unique.sorted.BL.bam.readcounts"
    shell:
        "samtools view -c {input} > {output}"

rule mergeBam:
    output: 
        "mergedBam/mergedSamples.bam"
    shell: 
        "cat filtered_sorted_BL_reads/*.bam > {output}"

rule bam2bed:
    input:
        "mergedBam/mergedSamples.bam"
    output:
        "mergedBam/mergedSamples.bed"
    shell:
        "bedtools bamtobed -i {input} > {output}"

rule peakcalling:
    input:
        "mergedBam/mergedSamples.bed"
    output:
        "peaks_macs/mergedSamples_macs_peaks.bed"
    shell:
        "macs2 callpeak -t {input} --keep-dup all -f BED -g hs -n mergedSamples --outdir peaks_macs -q 0.01"

       
