#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_PON.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create the sarek samplesheet for the WES normals used to build the PON, from
# the shared NAS samplesheet (the samples matching params$sample_regex).
#
# The fastqs live on the NAS and are NOT copied here: this script only writes
# the paths Snakemake will retrieve them to (params$nas_prefix + path on the
# NAS). The Sarek rule declares those files as nas() inputs, so each sample's
# fastqs are fetched when its job starts and deleted when it finishes.
#
# The samplesheet points at symlinks with sanitised names under
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
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input        <- snakemake@input[["samplesheet_nas"]]
    nas_prefix   <- snakemake@params[["nas_prefix"]]
    sample_regex <- snakemake@params[["sample_regex"]]
    links_dir    <- snakemake@params[["links_dir"]]
    output       <- snakemake@output[["samplesheet"]]
}else{
    input        <- '/trinity/home/r115502/NONC_Samplesheet/output/2026-10-21/Samplesheet_NAS.tsv'
    nas_prefix   <- '/trinity/home/r115502/PON/PON-snake/.snakemake/storage/nas/neuro-genomic-1.erasmusmc.nl'
    sample_regex <- 'Normal'
    links_dir    <- 'data/links'
    output       <- 'samplesheet.csv'
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

#-------------------------------------------------------------------------------
# 1.1 Read the NAS samplesheet and select this project's samples
#-------------------------------------------------------------------------------
Samplesheet_NAS <- read.delim(input)

Selection <- Samplesheet_NAS %>% filter(grepl(sample_regex, Sample))

if(nrow(Selection) == 0)
    stop("sample_regex '", sample_regex, "' matched no samples in ", input)

#-------------------------------------------------------------------------------
# 1.2 Build the sheet: one sample per patient (<patient>_normal1), one row per lane
#-------------------------------------------------------------------------------
rows <- list()

for(i in seq_len(nrow(Selection))){
    row     <- Selection[i, ]
    patient <- as.character(row$Sample)

    R1 <- split_list(row$R1_DNA)
    R2 <- split_list(row$R2_DNA)

    if(length(R1) == 0) next            # no sequencing data for this sample
    if(length(R1) != length(R2))
        stop("R1/R2 fastq counts differ for ", patient, ": ", length(R1), " vs ", length(R2))

    fastq_1 <- link_safe(retrieved(R1), patient)
    fastq_2 <- link_safe(retrieved(R2), patient)

    rows[[patient]] <- data.frame(
        patient = patient,
        sample  = paste0(patient, '_normal1'),   # the Snakefile's paths assume normal1
        status  = 1,                             # tumour-only mutect2 on each normal
        lane    = vapply(basename(fastq_1), function(b) strsplit(b, '_')[[1]][4],
                         character(1), USE.NAMES = FALSE),
        fastq_1 = fastq_1,
        fastq_2 = fastq_2)
}

if(length(rows) == 0)
    stop("none of the samples matching '", sample_regex, "' have fastq files")

Samplesheet <- bind_rows(rows) %>% arrange(patient, lane)

#-------------------------------------------------------------------------------
# 1.3 Write to file
#-------------------------------------------------------------------------------
dir.create(links_dir,       recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

write.table(Samplesheet, file = output, sep = ',', quote = F, row.names = F)
