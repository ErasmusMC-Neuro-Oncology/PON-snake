#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_MINT.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create samplesheet for WES analysis of MINT samples
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  30-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    data_dir <- snakemake@params[["data_dir"]]
    input_sample_overview <- snakemake@params[["sample_overview"]]
    output <- snakemake@output[[1]]
}else{
    data_dir <-  '~/mnt/BIGR_home/MINT/fastq/'
    input_sample_overview <- '../data/MINT_db.xlsx'
    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
# Read sample overview
sample_overview <- readxl::read_xlsx(input_sample_overview , skip = 2, sheet = 'DB_Iris')

# Fetch WES .fastq files
fastq_paths <- list.files(data_dir, pattern = '.fastq.gz', recursive = T, full.names = T)

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet
#-------------------------------------------------------------------------------
# Fetch fastq file IDs
fastq_files <- data.frame(fastq = fastq_paths) %>%
    mutate(GS_ID = purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][2]),
           read = ifelse(grepl('R1.fastq.gz', fastq),'fastq_1','fastq_2'),
           lane =  purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][4])) %>%
    tidyr::pivot_wider(names_from = read , values_from = fastq)

# Fetch samplenames
sample_overview <-
    sample_overview %>%
    mutate(patient = gsub('\\(res2\\)','_R2',paste0('MINT',gsub('^M','',`MINT subject nr`)))) %>%
    select(patient, GS_ID, Sentrix_ID) %>%
    group_by(patient) %>%
    mutate(sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("tumor", row_number())))) %>%
    ungroup() %>%
    mutate(status = 1)

#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
fastq_files %>% left_join(sample_overview) %>% mutate(status = 1) %>%  select(patient,sample,status,lane,fastq_1,fastq_2) %>% write.table(file=output, sep = ',', quote = F, row.names = F)
