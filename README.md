# PON-snake

## Overview

This pipeline creates the **normal reference resources** required for the tumor-only `WES-snake` workflow. It is designed for sequencing data from non-tumor control samples processed with the same assay as the tumor cohort.

The workflow combines:

* **nf-core/Sarek** for alignment and normal-only mutation calling with Mutect2
* **CopywriteR / QDNAseq** for generation of copy number reference profiles

The main goal is to generate high-quality reference files used for:

* artifact filtering during tumor-only SNV/indel calling
* normalization during tumor-only CNA calling

---

## Workflow Summary

### 1. Sample Preparation

The pipeline first creates a `samplesheet.csv` from input metadata files. This script has to be edited by the user.

### 2. Preprocessing and mutation calling with Sarek

For each normal sample:

* Reads are aligned to the reference genome
* Duplicate marking and recalibration are performed
* Variants are called using **Mutect2**
* BAM files are retained for CNA processing

### 3. Mutation Panel of Normals (PON)

All Mutect2 normal VCF files are merged and used to generate:

* A **Mutect2 Panel of Normals (PON)**

This file can be used in `WES-snake` to suppress recurrent technical artifacts and background calls.
The recommended number of normal samples for such a PON is ~40 samples.

### 4. Copy Number Reference Generation

Using normal BAM files:

* Off-target read depth profiling with **CopywriteR**
* Read count generation with **QDNAseq**
* Per-bin reference panels are created for configured bin sizes

These files are used in `WES-snake` to normalize tumor copy number profiles.

### 5. Reporting / Export

Reference files and intermediate QC objects are written to disk.

---

## Software Requirements

* **Snakemake 9.17.3**
* **nf-core/Sarek > 3.8.1**
* Conda / Mamba

---

## Running the Pipeline

### Option 1 — SLURM Cluster

```bash
snakemake --profile slurm/
```

### Option 2 — Local Workstation

```bash
snakemake --cores 1 --use-conda
```

---

## Main Output Files

### Final Mutation Panel of Normals

```text
PON/PON.vcf.gz
```

Primary final mutation reference file used by Mutect2 in tumor-only calling.

### Final Copy Number Reference Files

```text
PON/CNA_PON_{binsize}.Rds
```

Primary CNA normalization outputs generated for each configured bin size.

### Additional Outputs

#### Per-Sample QDNAseq Read Count Objects

```text
QDNAseq/{binsize}/{sample}/data/QDNAseq_Readcounts.Rds
```

Read count objects used to construct CNA reference panels.

#### Segmentation / QC Files

```text
QDNAseq/{binsize}/{sample}/data/QDNAseq_Segmented.Rds
QDNAseq/{binsize}/{sample}/data/QDNAseq_Segments.txt
QDNAseq/{binsize}/{sample}/plots/QDNAseq_segmented_profile.pdf
```

Useful for evaluating copy number neutrality of normal samples.

#### Intermediate PON Construction Files

```text
PON/normal_list.txt
PON/merged_normals.vcf.gz
```

Used during Mutect2 PON generation.

#### Sarek Variant Outputs

```text
sarek/{sample}/variant_calling/mutect2/{sample}_normal1/{sample}_normal1.mutect2.filtered.vcf.gz
```

Per-sample normal variant calls used as input for the final PON.

---

## Directory Structure

```text
output/
├── sarek/
├── QDNAseq/
├── copywriter/
└── PON/
```

---

## Notes on Reference Sample Selection

For best results, normal samples should:

* Use the same capture kit as tumor samples
* Use the same sequencing platform
* Use similar coverage depth
* Be processed with similar library preparation methods

Larger numbers of normal samples generally improve robustness of the reference files.

---

## Relationship to WES-snake

This pipeline generates the normal reference resources used by the tumor-only `WES-snake` workflow:

```text
PON/PON.vcf.gz
PON/CNA_PON_{binsize}.Rds
```

These files should be supplied in the `config.yaml` of `WES-snake`
