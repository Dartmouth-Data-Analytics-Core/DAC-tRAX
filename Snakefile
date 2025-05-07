#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# GDSC-tRAX v2 Pipeline (Claude Version)
# 
# Authors: Mike Martinez
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SET GLOBAL SCOPE PYTHON VARIABLES (EXECUTED BEFORE SNAKEMAKE)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import pandas as pd 

#----- Set config file
configfile: "config.yaml"

#----- Read in the sample data
samples_df = pd.read_table(config["sample_txt"], delimiter = ",").set_index("Sample_ID", drop = False)
sample_list = list(samples_df["Sample_ID"])
genome = config["genome"]


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# SNAKEMAKE RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#----- Final Rule
rule all:
    input:
        #----- Rule trimming outputs
        expand("trimming/{sample}.R1.trim.fastq.gz", sample = sample_list),
        expand("trimming/logs/{sample}.cutadapt.report", sample = sample_list),

        #----- Rule alignment outputs
        expand("alignment/{sample}.bam", sample = sample_list),
    output:
        "done.txt"
    conda: "rnaseq1"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        genome = config["genome"]
    shell:"""
    
        #----- Make dummy file
        touch done.txt
    
    """

#----- Rule to trim
rule trimming:
    output:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
        report = "trimming/logs/{sample}.cutadapt.report"
    conda: "cutadapt"
    resources: cpus="8", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        fastq_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        adapter_1 = config["adapter_1"],
        minlength = config["minlength"]
    shell: """
    
        #----- Run cutadapt
        cutadapt \
            -o {output.trim_1} \
            {params.fastq_1} \
            -m {params.minlength} \
            -a {params.adapter_1} \
            -j {resources.cpus} > {output.report} 2>&1
    """

#----- Rule to align samples to tRNA database (need db-trnatable.txt)
rule align:
    input:
        trim_1 = "trimming/{sample}.R1.trim.fastq.gz",
    output:
        bam = "alignment/{sample}.bam"
    conda: "trax_env"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb"
    params:
        sample = lambda wildcards: wildcards.sample,
        database = config["trna_db"],
        bt2_index = config["bt2_index"],
        maxMaps = config["maxMaps"],
        nPenalty = config["nPenalty"]
    shell: """
    
        #----- Run Bowtie2 alignment (removed the --reorder arg from original code to speed up...don't see the need with SE data)
        bowtie2 \
            -x {params.bt2_index} \
            -U {input.trim_1} \
            -k {params.maxMaps} \
            --very-sensitive \
            --np {params.nPenalty} \
            --ignore-qual \
            -p {resources.cpus} | \
            samtools view -bS -o {output.bam} -
    """

#rule filter_and_sort_alignment:
#    input:
#    output:
#    conda:
#    resources:
#    params:
#    shell: """ this is using choosemappings.py"""




