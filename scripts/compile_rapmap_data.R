library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(Biostrings)
library(ggplot2)
library(stringr)
library(dplyr)
library(parallel)


# gff_file <- "~/Box Sync/data/projects/carTcellAnalysis/annotation/carPlus.gtf"
# car_gtf <- import.gff3(gff_file)

gff_file <- "data/annotation/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)

get_coverage <- function(aln, gtf, seq, meta = NULL, smooth = 0, min.coverage = 0) {

    ## set boundaries of region to plot: approx. the range of all transcripts
    message("Calculating boundaries of region of interest...")
    zoom <- IRanges::reduce(gtf, min.gapwidth = 100e4)
    zoom <- zoom[which(seqnames(zoom) == seq)]
    if (smooth > 0) {
        ## calculated smoothed coverage
        message("Smoothing coverage in windows of size ", smooth, " bp...")
        bins <- seq(start(zoom), end(zoom), smooth)
        bins_gr <- GRanges(seqnames = as.vector(seqnames(zoom))[1], 
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
        lib_id <- str_extract(x, "lib[0-9]+")
        aln <- readGAlignments(x)
        
        cov_df <- get_coverage(aln, gtf, seq_name, smooth = smooth) %>% 
            mutate(lib_id = lib_id)
    }) %>% 
        bind_rows()
}

### Collect & store coverage data for all samples ###

bam_files <- list.files("data/results/rapmap/",
                        full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".bam$")]

# ...for CAR transcript
car_gtf <- xcripts_gtf %>% 
    subset(., elementMetadata(.)$gene_id == "CAR")
car_cov_dat <- build_coverage_df(as.list(bam_files[1]), car_gtf, "CAR-1", 10)

# egfr_xcripts <- xcripts_gtf %>% 
#     elementMetadata() %>% 
#     as.data.frame() %>% 
#     filter(str_detect(transcript_id, "EGFRt_")) %>% 
#     .$transcript_id
# 
# # ...for EGFR transcripts
# # (note: this is really slow - should consider using mclapply)
# egfr_cov_dat <- mclapply(as.list(egfr_xcripts), function(x) {
#     egfr_gtf <- xcripts_gtf %>% 
#         subset(., elementMetadata(.)$transcript_id == x)
#     seqname <- seqnames(egfr_gtf) %>% 
#         as.character()
#     
#     xcript_cov_dat <- build_coverage_df(as.list(bam_files), egfr_gtf, 
#                                         seqname, 10) %>% 
#         mutate(egfr_xcript = x)
# }) %>% 
#     bind_rows()
# 
# 
# ### Format CAR GTF data ###
# xcript_dat <- as.data.frame(xcripts_gtf)
# 
# car_dat <- xcript_dat %>% 
#     filter(seqnames == "CAR-1") %>% 
#     mutate(segment = factor(transcript_id, levels = transcript_id))
# 
# 
# ### Build EGFRt GTF-like dataframe indicating transmembrane region ###
# 
# # read in reference transcript FASTA
# seqs <- readDNAStringSet("data/sequence/hg38_CAR_transcripts.fa")
# 
# # extract sequences of transcripts in 'xcripts_gtf'
# chr_names <- xcripts_gtf %>% 
#     seqnames() %>% 
#     as.list() %>% 
#     unique()
# 
# xcript_matches <- lapply(chr_names, function(x) {
#     which(str_detect(names(seqs), x))
# }) %>% as.numeric()
# 
# xcript_seqs <- seqs[xcript_matches]
# 
# # identify start & end positions for EGFRt segment
# egfrt_region <- car_gtf %>% 
#     subset(., elementMetadata(.)$transcript_id == "EGFRt") %>% 
#     ranges()
# 
# # extract EGFRt subsequence
# egfrt_subseq <- subseq(xcript_seqs[1], 
#                        start = egfrt_region@start - 1,
#                        width = egfrt_region@width + 2)
# 
# # find and store overlapping regions in reference EGFR sequences
# egfrt_dat <- lapply(as.list(9:13), function(x) {
#     seq_match <- pairwiseAlignment(egfrt_subseq[[1]], xcript_seqs[[x]], 
#                                    type = "local")
#     match_start <- seq_match %>% 
#         subject() %>% 
#         start()
#     match_end <- seq_match %>% 
#         subject() %>% 
#         width() + match_start - 1
#     
#     egfrt_match <- data_frame(seqnames = chr_names[[x]],
#                               start = match_start,
#                               end = match_end,
#                               segment = "transmembrane")
# }) %>% 
#     bind_rows() %>% 
#     left_join(xcripts_gtf %>% 
#                   as.data.frame() %>% 
#                   select(seqnames, transcript_id)) %>% 
#     mutate(segment = factor(segment, levels = segment))
# 
# 
# 
# # save image
# save(car_cov_dat, egfr_cov_dat, car_dat, egfrt_dat, file = "data/sample_rapmap_data.RData")

