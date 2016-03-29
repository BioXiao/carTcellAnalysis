library(GenomicRanges)
library(rtracklayer)
library(GenomicAlignments)
library(ggplot2)
library(dplyr)
# library(ggbio)

gff_file <- "~/Box Sync/data/projects/carTcellAnalysis/annotation/carPlus.gtf"
car_gtf <- import.gff3(gff_file)

get_coverage_df <- function(bam_file, seq_name) {
    lib_id <- str_extract(bam_file, "lib[0-9]+")
    
    gal <- readGAlignments(bam_file)
    cov <- coverage(gal)
    
    seq_cov <- cov[[seq_name]]
    seq_len <- length(seq_cov)
    cov_df <- data_frame(pos = 1:seq_len,
                         cov = as.vector(window(seq_cov, 1, seq_len)))
    return(cov_df)
}

get_coverage <- function(aln, gtf, seq, meta = NULL, smooth = 0, min.coverage = 0) {

    ## set boundaries of region to plot: approx. the range of all transcripts
    message("Calculating boundaries of region of interest...")
    zoom <- reduce(gtf, min.gapwidth = 100e4)
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
    lapply(bam_file_list, function(x) {
        lib_id <- str_extract(x, "lib[0-9]+")
        aln <- readGAlignments(x)
        
        cov_df <- get_coverage(aln, gtf, seq_name, smooth = smooth) %>% 
            mutate(lib_id = lib_id)
    }) %>% 
        bind_rows()
}


plot_coverage <- function(cov_df, xlims) {
    cov_df %>% 
        ggplot(aes(x = pos, y = cov, colour = lib_id)) +
        geom_point(alpha = 0.5) +
        geom_smooth(se = FALSE) +
        scale_color_viridis(discrete = TRUE) +
        theme_bw()
}


bam_files <- list.files("~/Box Sync/data/projects/carTcellAnalysis/results/rapmap/",
                        full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".bam$")]

cov_df <- build_coverage_df(as.list(bam_files), car_gtf, "CAR-1", 10)
