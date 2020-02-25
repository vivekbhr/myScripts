infolder = "../dedup_bam"
samples = ['Bl6-BM-K36me3-CTCF05-i1', 'Bl6-BM-K36me3-CTCF05-i2']
clusters_tsv = ['Bl6-BM-K36me3-CTCF05-i1.txt', 'Bl6-BM-K36me3-CTCF05-i2.txt']
splitbam_path = '~/programs/myScripts/splitBam.py'
import csv
import os

def getClusters(tsv):
    cluster_dict = {}
    with open(tsv) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        for row in csv_reader:
            cluster_dict[row[0]] = row[1]
    clusters = set(x for x in cluster_dict.values())
    return(list(clusters))

clusters = getClusters(clusters_tsv[0])

rule all:
    input: expand("clusters_split/{sample}_{cluster}.bam", sample = samples, cluster = clusters),
           expand("clusters_split/{sample}_{cluster}.bw", sample = samples, cluster = clusters),
           expand("clusters_split/{cl}_merged.bam.bai", cl = clusters),
           expand("clusters_split/{cl}_merged.bw", cl = clusters)

## split bam by cluster
rule splitBAM:
    input:
        bam = infolder+"/{sample}.bam",
        bai = infolder+"/{sample}.bam.bai",
        clusters = '{sample}.txt'
    output:
        bams = temp(expand("clusters_split/{{sample}}_{cluster}.bam", cluster = clusters))
    wildcard_constraints:
        cluster="\d"
    params:
        outdir = "clusters_split/{sample}_",
        splitbam = splitbam_path
    threads: 1
    shell:
        "{params.splitbam} -i {input.bam} -c {input.clusters} -o {params.outdir}"

## index bam
rule index:
    input: "clusters_split/{sample}_{cluster}.bam"
    output: temp("clusters_split/{sample}_{cluster}.bam.bai")
    wildcard_constraints:
        cluster="\d"
    threads: 1
    shell:
        "samtools index {input}"

## make bigwig
rule bamcov:
    input:
        bam = "clusters_split/{sample}_{cluster}.bam",
        bai = "clusters_split/{sample}_{cluster}.bam.bai"
    output: "clusters_split/{sample}_{cluster}.bw"
    wildcard_constraints:
        cluster="\d"
    params:
        blacklist = "/hpc/hub_oudenaarden/vbhardwaj/annotations/blacklists/mm10-blacklist.v2.bed"
    threads: 20
    shell:
        "bamCoverage -p {threads} --skipNAs -bl {params.blacklist} -ignore X Y M "
        "--effectiveGenomeSize 2652783500 "
        "-b {input.bam} -o {output}"

## optionally, merge bw across samples
rule mergeBam:
    input:
         expand("clusters_split/{sample}_{{cl}}.bam", sample = samples)
    output: "clusters_split/{cl}_merged.bam"
    wildcard_constraints:
        cl="\d"
    threads: 10
    shell:
        "samtools merge -@ {threads} {output} {input}"

rule idx_mergeBam:
    input: "clusters_split/{cl}_merged.bam"
    output: "clusters_split/{cl}_merged.bam.bai"
    wildcard_constraints:
        cl="\d"
    threads: 1
    shell:
        "samtools index {input}"

rule bamcov_mergeBam:
    input:
        bam = "clusters_split/{cl}_merged.bam",
        bai = "clusters_split/{cl}_merged.bam.bai"
    output: "clusters_split/{cl}_merged.bw"
    wildcard_constraints:
        cl="\d"
    params:
        blacklist = "/hpc/hub_oudenaarden/vbhardwaj/annotations/blacklists/mm10-blacklist.v2.bed"
    threads: 20
    log: "{cl}_bamCoverage.log"
    shell:
        "bamCoverage -p {threads} --skipNAs -bl {params.blacklist} -ignore X Y M "
        "--effectiveGenomeSize 2652783500 "
        "--normalizeUsing RPGC "
        "-b {input.bam} -o {output} > {log} 2>&1"
