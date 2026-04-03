#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_PON.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create samplesheet for WES analysis of normal samples
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
#  03-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    data_dir_MINT <- snakemake@params[["data_dir_MINT"]]
    input_sample_overview_MINT <- snakemake@params[["sample_overview_MINT"]]
    output <- snakemake@output[[1]]
}else{
    data_dir_MINT <-   '~/mnt/BIGR_home/MINT/fastq_batch2/'
    input_sample_overview_MINT <- '~/mnt/BIGR_home/MINT/data/1kuvre_fastq_list.csv'
    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
# Read sample overview
sample_overview_MINT <- read.delim(input_sample_overview_MINT, sep = ',')

# Fetch WES .fastq files
fastq_paths_batch_MINT <- list.files(data_dir_MINT, pattern = '.fastq.gz', recursive = T, full.names = T)

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet\
#-------------------------------------------------------------------------------
fastq_files_MINT <- data.frame(fastq = fastq_paths_batch_MINT) %>%
    mutate(GS_ID = purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][2]),
           read = ifelse(grepl('R1.fastq.gz', fastq),'fastq_1','fastq_2'),
           lane =  purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][4])) %>%
    tidyr::pivot_wider(names_from = read , values_from = fastq) %>%
    mutate(basenameR1 = basename(fastq_1))

sample_overview_MINT <- sample_overview_MINT %>%
    mutate(patient = gsub('_M','',purrr::map_chr(RGSM,~strsplit(.x,'-')[[1]][4])),
           sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("normal", row_number()))),
           basenameR1 = basename(Read1File)) %>%
    group_by(patient) %>%
    mutate(sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("normal", row_number())))) %>%
    ungroup() %>%
    filter(grepl('Normal',patient)) %>%
    mutate(status = 1)

print(sample_overview_MINT)
print(fastq_files_MINT)

sample_overview_MINT %>% left_join(fastq_files_MINT) %>% print()

#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
sample_overview_MINT %>% left_join(fastq_files_MINT) %>%  select(patient,sample,status,lane,fastq_1,fastq_2) %>% write.table(file=output, sep = ',', quote = F, row.names = F)
