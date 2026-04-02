configfile: "config.yaml"
from datetime import datetime
import pandas as pd
import os
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare variables and wildcards
output_dir = config["all"]["output_dir"]

# Fetch Sample wildcards
Samples = pd.read_csv('samplesheet.csv')['patient'].unique() if os.path.isfile('samplesheet.csv') else []

#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        'samplesheet.csv',
        expand(output_dir + 'QDNAseq/{binsize}/{sample}/plots/QDNAseq_segmented_profile.pdf', sample = Samples, binsize = config['CopyWriteR']['binsizes'])

#++++++++++++++++++++++++++++++++++++++++++++ 0 CREATE SAMPLESHEET ++++++++++++++++++++++++++++++++++++++++++++++++
rule Create_Samplesheet:
    params:
        data_dir = config['all']['data_dir'],
        sample_overview = '../data/MINT_db.xlsx',
    output:
        'samplesheet.csv'
    conda:
        'envs/R.yaml'
    script:
        'scripts/Create_Samplesheet_PON.R'

        
#+++++++++++++++++++++++++++++++++++++++++ 1 RUN SAREK VARIANT CALLING +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Download Sarek
rule Download_Sarek:
    output:
        ".nf-core-sarek/"
    params:
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/"
    conda:
        "envs/nextflow.yaml"
    shell:
        """
	export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}        
	nf-core pipelines download --outdir {output} --container-system singularity --compress none -r dev sarek
        """

# 1.2 Run Sarek
rule Sarek:
    input:
        samplesheet = 'samplesheet.csv',
        sarek = ".nf-core-sarek/"
    output:
        samplesheet = output_dir + 'sarek/{sample}/csv/samplesheet.csv',
        vcf = 'output/sarek/{sample}/variant_calling/mutect2/{sample}_normal1/{sample}_normal1.mutect2.filtered.vcf.gz',
        mapped = temp(directory(output_dir + "sarek/{sample}/preprocessing/mapped/")),
        recal = temp(directory(output_dir + "sarek/{sample}/preprocessing/recalibrated/")),
        md_cram = temp(output_dir + "sarek/{sample}/preprocessing/markduplicates/{sample}_normal1/{sample}_normal1.md.cram"),
        md_bam = temp(output_dir + "sarek/{sample}/preprocessing/markduplicates/{sample}_normal1/{sample}_normal1.md.bam"),
    threads: 8
    resources:
        mem_mb=50000,
        gpu=0,
        runtime='30h'
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/sarek/{sample}/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H%M%S")+".log"
    params:
        version = 'dev',
        genome = 'GATK.GRCh38',
        profile = "singularity",
        tools = "mutect2",
        Mutect2_params = 'params/mutect2_params.json',
        targets = config['sarek']['targetregions'],
        intervals = config['sarek']['interval_padding'],
        singularity_dir = f"{config['sarek']['workdir']}/singularity/cache/",
        workdir=lambda wildcards: f"{config['sarek']['workdir']}/{wildcards.sample}",
        outdir=lambda wildcards: f"{output_dir}/sarek/{wildcards.sample}",
    shell:
        """
        export NXF_WORK={params.workdir}
        export NXF_SINGULARITY_CACHEDIR={params.singularity_dir}

        # Subset samplesheet
        awk -F',' '$1=="sample" || $1=="{wildcards.sample}"' {input.samplesheet} > {output.samplesheet}

        nextflow -log {log} run .nf-core-sarek/{params.version}/ \
           -profile {params.profile} \
           -work-dir {params.workdir} \
           -c {params.Mutect2_params} \
              --input {output.samplesheet} \
              --outdir {params.outdir} \
              --genome {params.genome} \
              --tools {params.tools} \
              --intervals {params.targets} \
              --interval_padding {params.intervals} \
              --max_memory '{resources.mem_mb} MB' \
              --save_mapped  \
              --wes

        # Save alignment as .bam (to be fixed with --save-output-as-bam in new sarek release)
        samtools view -b -o {output.md_bam} {output.md_cram}
        
        # Clean cache and intermediate files upon completion but keep on failure
        status=$?
        if [ $status -eq 0 ]; then
        rm -rf {params.workdir}
        else
        echo "Sarek failed"
        fi
        
        """


#-------------------------------------------------------------------------------------------------------------        
# 1.2 Create Mutect2 PON
rule Create_PON:
     input:
        expand('output/sarek/{sample}/variant_calling/mutect2/{sample}_normal1/{sample}_normal1.mutect2.filtered.vcf.gz', sample = Samples)
    output:
        merged_vcf = output_dir + 'PON/merged_normals.vcf.gz'
        PON = output_dir + 'PON/PON.vcf.gz'
    threads: 8
    resources:
        mem_mb=50000,
        gpu=0,
    conda:
        "envs/nextflow.yaml"
    params:
        ref = config['all']['ref'],
    shell:
        """
        # merge vcfs
        bcftools merge {input} --merge all -o {output.merged_vcf}
        
        gatk CreateSomaticPanelOfNormals \
        -R {params.ref} \
        -V {output.merged_vcf} \
        -O {output}
        """


#+++++++++++++++++++++++++++++++++++++++++ 2 PERFORM CNA ANALYSIS +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Run CopywriteR, QDNAseq, ACE, CNH, calculate stats and export results
rule CNA_analysis:
    input:
        bam= output_dir + "sarek/{sample}/preprocessing/markduplicates/{sample}_normal1/{sample}_normal1.md.bam"
    output:
        sample_dir = temp(directory(output_dir + "copywriter/{binsize}/{sample}/")),
        QDNAseq = output_dir + 'QDNAseq/{binsize}/{sample}/data/QDNAseq_Segmented.Rds',
        Profile = output_dir + 'QDNAseq/{binsize}/{sample}/plots/QDNAseq_segmented_profile.pdf',
        Segments = output_dir + 'QDNAseq/{binsize}/{sample}/data/QDNAseq_Segments.txt',
        Segments_igv = temp(output_dir + 'QDNAseq/{binsize}/{sample}/data/QDNAseq_Segments.igv'),
        Called = output_dir + 'QDNAseq/{binsize}/{sample}/data/QDNAseq_calls.txt',
        CNA_stats = output_dir + 'QDNAseq/{binsize}/{sample}/data/CNA_stats.txt',
        ACE_results = output_dir + 'ACE/{binsize}/{sample}/ACE_fits.txt',
        ACE_matrix = output_dir + 'ACE/{binsize}/{sample}/ACE_matrixplot.pdf',
        CNH_results = output_dir + 'CNH/{binsize}/{sample}/CNH_results.txt',
        CNH_plot = output_dir + 'CNH/{binsize}/{sample}/CNH_plot.pdf',
        CNH_error_plot = output_dir + 'CNH/{binsize}/{sample}/CNH_errorplot.pdf',
    resources:
        mem_mb=10000,
        gpu=0,
        runtime='30h'
    params:
        genome = 'hg38',
        cores = config['CopyWriteR']['cores'],
        cytobands = config['CopyWriteR']['cytobands'],
        ACE_purity_penalty = config['ACE']['penalty'],
        ACE_ploidy_penalty = config['ACE']['penploidy'],
        CNH_path = config['CNH']['path']
    conda:
        "envs/CNA.yaml"
    script:
        'scripts/CNA_analysis.R'
