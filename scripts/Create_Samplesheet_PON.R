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
    data_dir <- snakemake@params[["data_dir"]]
    output <- snakemake@output[[1]]
}else{
    data_dir <-   '~/mnt/BIGR_home/PON/PON_Heidelberg_panel/fastq/'
    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
# Fetch WES .fastq files
fastq_paths <- list.files(data_dir, pattern = '.fastq.gz', recursive = T, full.names = T)


#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet\
#-------------------------------------------------------------------------------
fastq_files <- data.frame(fastq = fastq_paths) %>%
    mutate(patient = purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][1]),
           sample = paste0(patient,'_normal1'),
           read = ifelse(grepl('R1_001.fastq.gz', fastq),'fastq_1','fastq_2'),
           lane =  purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][3])) %>%
    tidyr::pivot_wider(names_from = read , values_from = fastq) %>%
    mutate(status=1)


#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
fastq_files %>%  select(patient,sample,status,lane,fastq_1,fastq_2) %>% write.table(file=output, sep = ',', quote = F, row.names = F)
