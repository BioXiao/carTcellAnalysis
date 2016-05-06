
# load packages -----------------------------------------------------------

library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(Biostrings)
library(stringr)
library(dplyr)
library(parallel)

# load data ---------------------------------------------------------------

gff_file <- "data/annotation/_old/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)
xcript_dat <- as.data.frame(xcripts_gtf) %>% 
    mutate(seqnames = as.character(seqnames))

# load functions ----------------------------------------------------------

source("scripts/data_cleaning_functions.R")

# define functions --------------------------------------------------------

get_coverage <- function(aln, gtf, seq, smooth = 0, min.coverage = 0) {

    ## set boundaries of region to plot: approx. the range of all transcripts
    message("Calculating boundaries of region of interest...")
    zoom <- IRanges::reduce(gtf, min.gapwidth = 100e4)
    zoom <- zoom[which(seqnames(zoom) == seq)]
    if (smooth > 0) {
        ## calculating smoothed coverage
        message("Smoothing coverage in windows of size ", smooth, " bp...")
        bins <- seq(start(zoom), end(zoom), smooth)
        bins_gr <- GRanges(seqnames = as.vector(seqnames(zoom)), 
                           ranges = IRanges(start = bins [-length(bins)], end = bins[-1]-1))

        cov_df <- transform(as.data.frame(bins_gr), 
                            score = countOverlaps(bins_gr, aln)) %>% 
            mutate(midpoint = (start + end) / 2) %>% 
            select(pos = midpoint, cov = score)
    }
    else {
        ## calculate raw coverage
        message("Estimating read coverage...")
        cvg <- coverage(aln)[[seq]]
        
        seq_len <- length(cvg)
        cov_df <- data_frame(pos = 1:seq_len,
                             cov = as.vector(window(cvg, 1, seq_len)))
    }
    
    return(cov_df)
}

build_coverage_df <- function(bam_file_list, gtf, seq_name, smooth = 0) {
    mclapply(bam_file_list, function(x) {
        x_lib_id <- str_extract(x, "lib[0-9]+")
        aln <- readGAlignments(x)
        
        cov_df <- get_coverage(aln, gtf, seq_name, smooth = smooth) %>% 
            mutate(lib_id = x_lib_id)
    }) %>% 
        bind_rows()
}

# specify file paths ------------------------------------------------------

bam_files <- list.files("data/results/_old/rapmap",
                        full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".bam$")] %>% 
    .[str_detect(., "lib8845")]

# collect & format CAR coverage data for all samples ----------------------

car_gtf <- xcripts_gtf %>% 
    subset(., elementMetadata(.)[["gene_id"]] == "CAR")
car_cov_dat <- build_coverage_df(as.list(bam_files), car_gtf, "CAR-1", 10)

# add segment to CAR coverage data ----------------------------------------

# function to append CAR segment
get_segment <- Vectorize(function(pos, seq) {
    pos <- ceiling(pos)
    xcript_dat %>% 
        filter(seqnames == seq,
               (start - 1) <= pos & (end + 1) >= pos) %>% 
        select(transcript_id) %>% 
        as.character()
})

lib_cov <- car_cov_dat %>% 
    filter(lib_id == car_cov_dat[["lib_id"]][1]) %>% 
    mutate(segment = get_segment(pos, "CAR-1")) %>% 
    select(pos, segment)

car_cov_dat <- car_cov_dat %>% 
    left_join(lib_cov, by = c("pos" = "pos")) %>% 
    filter(segment != "character(0)")

# collect & format EGFR coverage for all samples --------------------------

egfr_xcripts <- xcripts_gtf %>% 
    elementMetadata() %>% 
    as.data.frame() %>% 
    filter(str_detect(transcript_id, "EGFRt_")) %>% 
    .[["transcript_id"]]

egfr_cov_dat <- mclapply(as.list(egfr_xcripts), function(x) {
    egfr_gtf <- xcripts_gtf %>% 
        subset(., elementMetadata(.)[["transcript_id"]] == x)
    seqname <- seqnames(egfr_gtf) %>% 
        as.character()
    
    xcript_cov_dat <- build_coverage_df(as.list(bam_files), egfr_gtf, 
                                        seqname, 10) %>% 
        mutate(egfr_xcript = x)
}) %>% 
    bind_rows()

# format CAR GTF data -----------------------------------------------------

car_dat <- xcript_dat %>% 
    filter(seqnames == "CAR-1") %>% 
    mutate(segment = factor(transcript_id, levels = transcript_id))

# format EGFR GTF data ----------------------------------------------------

# read in reference transcript FASTA
seqs <- readDNAStringSet("data/sequence/_old/hg38_CAR_transcripts.fa")

# extract sequences of CAR & EGFR transcripts in 'xcripts_gtf'
chr_names <- xcripts_gtf %>% 
    subset(., elementMetadata(.)[["gene_id"]] %in% c("CAR", "EGFR")) %>% 
    seqnames() %>% 
    as.list() %>% 
    unique()

xcript_matches <- lapply(chr_names, function(x) {
    which(str_detect(names(seqs), x))
}) %>% as.numeric()

xcript_seqs <- seqs[xcript_matches]

# identify start & end positions for EGFRt segment
egfrt_region <- car_gtf %>% 
    subset(., elementMetadata(.)$transcript_id == "EGFRt") %>% 
    ranges()

# extract EGFRt subsequence
egfrt_subseq <- subseq(xcript_seqs[1], 
                       start = start(egfrt_region) - 1,
                       width = width(egfrt_region) + 2)

# find and store overlapping regions in reference EGFR sequences
egfrt_dat <- lapply(as.list(2:length(chr_names)), function(x) {
    seq_match <- pairwiseAlignment(egfrt_subseq[[1]], xcript_seqs[[x]], 
                                   type = "local")
    match_start <- seq_match %>% 
        subject() %>% 
        start()
    match_end <- seq_match %>% 
        subject() %>% 
        width() + match_start - 1
    
    egfrt_match <- data_frame(seqnames = chr_names[[x]],
                              start = match_start,
                              end = match_end,
                              segment = "EGFRt")
}) %>% 
    bind_rows() %>% 
    mutate(seqnames = as.character(seqnames)) %>% 
    left_join(xcript_dat %>% 
                  select(seqnames, transcript_id)) %>% 
    mutate(segment = factor(segment, levels = c("EGFRt", "other"))) %>% 
    dplyr::rename(egfr_xcript = transcript_id)

# save compiled files -----------------------------------------------------

write_csv(car_cov_dat, "data/clean/all_compiled_car_coverage_data.csv")
write_csv(car_dat, "data/clean/car_gtf_data.csv")

write_csv(egfr_cov_dat, "data/clean/all_compiled_egfr_coverage_data.csv")
write_csv(egfrt_dat, "data/clean/egfr_gtf_data.csv")

