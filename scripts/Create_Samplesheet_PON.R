#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_PON.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create the sarek samplesheets for the WES normals used to build the PON, from
# the shared NAS samplesheet (the samples matching params$sample_regex).
#
# ONE SHEET PER SAMPLE, because the two input kinds need different columns and
# each sample is its own sarek run:
#   R1_DNA = *.fastq.gz -> patient,sample,status,lane,fastq_1,fastq_2  (--step mapping)
#   R1_DNA = *.bam      -> patient,sample,status,bam,bai               (--step markduplicates)
#
# The files live on the NAS and are NOT copied here: this script only writes the
# paths Snakemake will retrieve them to (params$nas_prefix + path on the NAS).
# The Sarek rule declares them as nas() inputs, so each sample's files are
# fetched when its job starts and deleted when it finishes.
#
# The sheets point at symlinks with sanitised names under
# <data_dir>/links/<patient>/, because NAS paths may contain spaces and
# brackets, which nf-core samplesheet validation rejects. A symlink may dangle
# until the file is fetched; it resolves by the time sarek reads it.
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: R
#
# History:
#  03-04-2026: File creation
#  30-10-2026: fastqs are read from the NAS, samples from the NAS samplesheet
#  02-11-2026: Accept pre-aligned BAM input; one sheet per sample
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input           <- snakemake@input[["samplesheet_nas"]]
    nas_prefix      <- snakemake@params[["nas_prefix"]]
    sample_regex    <- snakemake@params[["sample_regex"]]
    links_dir       <- snakemake@params[["links_dir"]]
    samplesheet_dir <- snakemake@params[["samplesheet_dir"]]
}else{
    input           <- '/trinity/home/r115502/NONC_Samplesheet/output/2026-10-21/Samplesheet_NAS.tsv'
    nas_prefix      <- '/trinity/home/r115502/PON/PON-snake/.snakemake/storage/nas/neuro-genomic-1.erasmusmc.nl'
    sample_regex    <- 'Normal'
    links_dir       <- 'data/links'
    samplesheet_dir <- 'samplesheets'
}

#-------------------------------------------------------------------------------
# 0.3 Helpers
#-------------------------------------------------------------------------------
# ';'-separated NAS paths -> character vector ('NA'/empty -> character(0))
split_list <- function(x){
    if(length(x) == 0 || is.na(x)) return(character(0))
    v <- trimws(strsplit(as.character(x), ';')[[1]])
    v[nzchar(v) & v != 'NA']
}

# where Snakemake retrieves a NAS file to
retrieved <- function(relpath)
    file.path(nas_prefix, sub('^/', '', trimws(relpath)))

# symlink with a shell-safe name, in this patient's link directory
link_safe <- function(paths, patient){
    if(length(paths) == 0) return(character(0))
    pdir <- file.path(links_dir, patient)
    dir.create(pdir, recursive = TRUE, showWarnings = FALSE)
    names_ok <- gsub('[^A-Za-z0-9._-]', '_', basename(paths))
    if(anyDuplicated(names_ok))
        stop("sanitised file names collide for ", patient, ": ",
             paste(names_ok[duplicated(names_ok)], collapse = ', '))
    dest <- file.path(normalizePath(pdir), names_ok)
    for(i in seq_along(paths)){
        unlink(dest[i])
        # the target may not exist yet (it is fetched when the Sarek job runs)
        if(!file.symlink(normalizePath(paths[i], mustWork = FALSE), dest[i]))
            stop("could not link ", paths[i], " -> ", dest[i])
    }
    dest
}

is_bam <- function(x) length(x) > 0 && grepl('\\.bam$', x[1], ignore.case = TRUE)

#-------------------------------------------------------------------------------
# 1.1 Read the NAS samplesheet and select this project's samples
#-------------------------------------------------------------------------------
Samplesheet_NAS <- read.delim(input)

Selection <- Samplesheet_NAS %>% filter(grepl(sample_regex, Sample))

if(nrow(Selection) == 0)
    stop("sample_regex '", sample_regex, "' matched no samples in ", input)

#-------------------------------------------------------------------------------
# 1.2 Write one sheet per sample
#-------------------------------------------------------------------------------
dir.create(links_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(samplesheet_dir, recursive = TRUE, showWarnings = FALSE)

n_written <- 0

for(i in seq_len(nrow(Selection))){
    row     <- Selection[i, ]
    patient <- as.character(row$Sample)

    R1 <- split_list(row$R1_DNA)
    R2 <- split_list(row$R2_DNA)

    if(length(R1) == 0) next            # no sequencing data for this sample

    if(is_bam(R1)){
        # --- pre-aligned BAM: sarek --step markduplicates ---
        if(length(R1) != 1)
            stop("expected exactly 1 BAM for ", patient, ", found ", length(R1))
        if(length(R2) > 0)
            stop("R2_DNA must be empty for BAM input (", patient, ")")
        Sheet <- data.frame(
            patient = patient,
            sample  = paste0(patient, '_normal1'),
            status  = 1,
            bam     = link_safe(retrieved(R1), patient),
            bai     = link_safe(retrieved(paste0(R1, '.bai')), patient))
    }else{
        # --- fastq: sarek --step mapping ---
        if(length(R1) != length(R2))
            stop("R1/R2 fastq counts differ for ", patient, ": ",
                 length(R1), " vs ", length(R2))
        fastq_1 <- link_safe(retrieved(R1), patient)
        fastq_2 <- link_safe(retrieved(R2), patient)
        Sheet <- data.frame(
            patient = patient,
            sample  = paste0(patient, '_normal1'),
            status  = 1,
            lane    = vapply(basename(fastq_1), function(b) strsplit(b, '_')[[1]][4],
                             character(1), USE.NAMES = FALSE),
            fastq_1 = fastq_1,
            fastq_2 = fastq_2)
    }

    write.table(Sheet, file = file.path(samplesheet_dir, paste0(patient, '.csv')),
                sep = ',', quote = F, row.names = F)
    n_written <- n_written + 1
}

if(n_written == 0)
    stop("none of the samples matching '", sample_regex, "' have sequencing data")
