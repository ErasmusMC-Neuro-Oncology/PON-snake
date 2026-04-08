#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CNA_analysis.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run CopywriteR, Normalize with QDNAseq and export results
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
# 1) Add PON corrections
#
# History:
#  13-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
if(!'QDNAseq.hg38' %in% installed.packages()){devtools::install_github("asntech/QDNAseq.hg38@main")}
suppressMessages(library(dplyr))
suppressMessages(library(CopywriteR))
suppressMessages(library(QDNAseq))
suppressMessages(library(QDNAseq.hg38))
suppressMessages(library(Biobase))
suppressMessages(library(GenomicRanges))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_bam <- snakemake@input[["bam"]]
    sample <- snakemake@wildcards[["sample"]]
    genome <- snakemake@params[["genome"]]
    binsize <- snakemake@wildcards[["binsize"]]
    cores <- snakemake@params[["cores"]]
    sample_dir <- snakemake@output[["sample_dir"]]
    QDNAseq_output <- snakemake@output[["QDNAseq"]]
    Profile_output <- snakemake@output[["Profile"]]  
    Segments_output <- snakemake@output[["Segments"]]
    
}else{
    input_bam <- '../output/sarek/Normal_brain_2/preprocessing/markduplicates/Normal_brain_2_normal1/Normal_brain_2_normal1.md.bam'
    sample <- 'Normal_brain_2_normal1'
    sample_dir <- '../output/copywriter/1000kbp/Normal_brain_2/'
    QDNAseq_output <- '../output/QDNAseq/1000kbp/Normal_brain_2/data/QDNAseq_Segments.Rds'
    Segments_output <- '../output/QDNAseq/1000kbp/Normal_brain_2/data/QDNAseq_Segments.txt'
    Profile_output <- '../output/QDNAseq/1000kbp/Normal_brain_2/plots/QDNAseq_segmented_profile.pdf'
    genome <- 'hg38'
    binsize <- '1000kbp'
    cores <- 1
    
}

#-------------------------------------------------------------------------------
# 1.1 Define CopyWritR parameters
#-------------------------------------------------------------------------------
binsize <- format(as.integer(gsub('kbp','',binsize))*1000, scientific = F)

# Create annotation files
if(!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
preCopywriteR(output.folder = sample_dir,
              bin.size = as.integer(binsize),
              ref.genome = genome,
              prefix = "chr")

# get number of kb bins
kbbin <- substring(binsize,1,nchar(binsize)-3)
# Load annotation files
load(file = file.path(sample_dir, paste0(genome,"_",kbbin,"kb_chr"), "blacklist.rda"))
# set number of workers
bp.param <- SnowParam(workers = cores, type = "SOCK")
# Create sample df
sample.control <- data.frame(samples = input_bam,controls=input_bam)

#-------------------------------------------------------------------------------
# 2.1 Run CopyWritR
#-------------------------------------------------------------------------------
# Run CopyWriteR
CopywriteR(sample.control = sample.control,
           destination.folder = sample_dir,
           reference.folder = file.path(sample_dir, paste0(genome,"_",kbbin,"kb_chr")),
           bp.param=bp.param)


#-------------------------------------------------------------------------------
# 3.1 Parse Copywriter output
#-------------------------------------------------------------------------------
read_counts <- read.delim(paste0(sample_dir,'/CNAprofiles/read_counts.txt'))
#-------------------------------------------------------------------------------
# 3.2 Prepare QDNAseq objects
#-------------------------------------------------------------------------------
#---------- create bins file  ----------
kbbin <- substring(binsize,1,nchar(binsize)-3)
load(paste0(sample_dir,'/',genome,"_",kbbin,"kb_chr/GC_mappability.rda"))

# create dataframe containing fdata fields
fData_all <-
    cbind(as.data.frame(seqnames(GC.mappa.grange)),
          as.data.frame(ranges(GC.mappa.grange)),
          as.data.frame((mcols(GC.mappa.grange))))
# create features as row names
rownames(fData_all) <- paste0(fData_all$value,":",fData_all$start,"-",fData_all$end)

bins <- getBinAnnotations(as.integer(kbbin), genome="hg38")

pData_bins <- pData(bins)
features_bins <- paste0('chr',rownames(pData_bins))
# remove weird features from bin data
pData_bins <- pData_bins[!rownames(pData_bins) %in% features_bins[!features_bins %in% rownames(fData_all)], ]
bins  <-pData_bins
#--------------------------------------------
#---------- create copynumber file ----------
colnames(read_counts)[6] <- 'read_counts'
read_counts$samples <- sample

counts  <-
    read_counts %>% 
    dplyr::select(Feature,read_counts,samples) %>%
    tidyr::spread(samples,read_counts,fill=0)

counts <- counts[match(paste0('chr',rownames(bins)), counts$Feature),] %>%
    dplyr::group_by(Feature) %>% dplyr::filter(dplyr::row_number() == 1 ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(Feature)) %>% 
    tibble::column_to_rownames(.,"Feature") %>%
    #tibble::add_column(bins = rownames(.), .before = colnames(.)[1]) %>%
    as.matrix()

#-------------------------------------------
#---------- create phenodata file ----------
phenodata <- data.frame(name = sample,
           total.reads = colSums(counts),
           used.reads = colSums(counts),
           row.names = sample,
           expected.variance = 0)

#--------------------------------------------------
#--------- create QDNAseqReadCounts object --------
rownames(bins) <- paste0('chr',rownames(bins))
features <- intersect(rownames(bins),rownames(counts))

QDNAseqCopyNumbers <- new("QDNAseqReadCounts",bins=bins[features,],counts=as.matrix(counts[features,]),phenodata=phenodata)

#-------------------------------------------------------------------------------
# 4.1 Perform QDNAseq normalizations
#-------------------------------------------------------------------------------
corrected <- applyFilters(QDNAseqCopyNumbers, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE , chromosomes=c('chrY','chrX','X','Y')) %>%
    estimateCorrection() %>%
    correctBins() %>%
    normalizeBins() %>%
    smoothOutlierBins() %>%
    segmentBins() %>%
    normalizeSegmentedBins()

#-------------------------------------------------------------------------------
# 4.2 Plot QDNAseq profile and callBins
#-------------------------------------------------------------------------------
pdf(Profile_output, width = 6 , height = 5)
plot(corrected)
dev.off()

#-------------------------------------------------------------------------------
# 4.3 Retrieve segments
#-------------------------------------------------------------------------------
# Fetch segments and calculate their values
Segments <- fData(corrected) %>%
    mutate(seg.mean =assayData(corrected)$segmented[,1]) %>%
    filter(!is.na(seg.mean)) %>%
    arrange(chromosome, start) %>%
  group_by(chromosome) %>%
  mutate(group = cumsum(seg.mean != lag(seg.mean, default = dplyr::first(seg.mean)))) %>%
  group_by(chromosome, group) %>%
  summarise(
    loc.start = min(start),
    loc.end   = max(end),
    num.mark  = n(),
    seg.mean  = dplyr::first(seg.mean),
    .groups = "drop"
  ) %>%
    mutate(
        chrom_order = case_when(
            chromosome %in% c("X","x") ~ 23,
            chromosome %in% c("Y","y") ~ 24,
            TRUE ~ as.numeric(chromosome)),
        ID = sampleNames(corrected),
        chrom = paste0('chr',chromosome)) %>%
  arrange(chrom_order, loc.start) %>%
    select(ID, chrom, loc.start, loc.end, num.mark, seg.mean) 



#-------------------------------------------------------------------------------
# 5.1 Write to file
#-------------------------------------------------------------------------------
# Save QDNAseq object
saveRDS(corrected,QDNAseq_output)

# Save segments
write.table(Segments, Segments_output, sep = '\t',quote = F, row.names = F)
