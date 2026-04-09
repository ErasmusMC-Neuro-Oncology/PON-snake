#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_CNA_PON.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine normal sampels
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
#  07-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(QDNAseq))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_QDNAseq <- snakemake@input[["QDNAseq"]]
    output_PON <- snakemake@output[["PON"]]
   
}else{
    input_QDNAseq <- Sys.glob('../output/QDNAseq/100kbp/*/data/QDNAseq_Segmented.Rds')
    output_PON <- '../output/PON/CNA_PON_100kbp.Rds'
    
}

#-------------------------------------------------------------------------------
# 1.1 Read and combine QDNAseq objects
#-------------------------------------------------------------------------------
for(i in seq_along(input_QDNAseq)){
    if(i == 1){
        CNA_data <- readRDS(input_QDNAseq[i])
    }else{
        CNA_data <- Biobase::combine(CNA_data, readRDS(input_QDNAseq[i]))
    }
}


#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
# Save QDNAseq object
saveRDS(PON,output_PON)
