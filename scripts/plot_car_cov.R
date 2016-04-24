
# load libs ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(ggthemes)
library(viridis)


# load data ---------------------------------------------------------------

load("data/sample_metrics_data.RData")
load("data/sample_rapmap_data.RData")
load("data/sample_salmon_data.RData")
# load("data/sample_tcr_data.RData")

gff_file <- "data/annotation/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)



# define functions --------------------------------------------------------

# create/format data frame for plotting coverage
prep_cov_dat <- function(lib_dat, cov_dat, metric_dat) {
    lib_dat %>% 
        clean_dup_names() %>% 
        simplify_lib_id() %>%
        clean_donor_ids() %>% 
        select(lib_id, donor_id, timepoint) %>% 
        relabel_timepoints() %>% 
        left_join(cov_dat, by = c("lib_id" = "lib_id")) %>% 
        left_join(metric_dat %>% 
                      clean_dup_names() %>% 
                      simplify_lib_id() %>% 
                      select(lib_id, fastq_total_reads, 
                             median_cv_coverage, mapped_reads_w_dups),
                  by = c("lib_id" = "lib_id")) %>% 
        left_join(map_rate_dat, by = c("lib_id" = "lib_id")) %>% 
        mutate(mapped_reads = fastq_total_reads * map_rate,
               norm_cov = cov / mapped_reads) %>% 
        left_join(salmon_quant_dat %>% 
                      select(lib_id, Name, TPM) %>% 
                      spread(Name, TPM) %>% 
                      dplyr::rename(CAR = `CAR-1`) %>% 
                      relabel_transcripts(xcript_dat))
}

# crazy function to build plots
plot_coverage <- function(formatted_cov_dat, gtf_dat, 
                          split = TRUE, fits_only = FALSE, 
                          color_by = "lib_num", fade_by = "log2(CAR + 1)",
                          hide_legend = FALSE) {
    
    # add lib numbers for coloring
    formatted_cov_dat <- formatted_cov_dat %>% 
        group_by(donor_id, timepoint) %>% 
        mutate(lib_num = dense_rank(lib_id) / n_distinct(lib_id),
               lib_num = as.character(lib_num)) %>% 
        ungroup()
    
    # determine plot height
    height <- log2(max(formatted_cov_dat$cov, na.rm = TRUE) + 1)
    
    # create color scale
    if (!fits_only) {
        gradient_stop <- myCbPal[3]
    } else {
        gradient_stop <- myCbPal[2]
    }
    num_libs <- n_distinct(formatted_cov_dat$lib_id)
    cc <- seq_gradient_pal(myCbPal[3], gradient_stop)(seq(0, 1, length.out = num_libs))
    
    # build plot
    p_cov <- ggplot() +
        geom_rect(data = gtf_dat, 
                  aes(xmin = start, xmax = end, ymin = 0, ymax = height, 
                      fill = segment),
                  alpha = 0.5, colour = "gray")
    if (!fits_only) {
        p_cov <- p_cov +
            geom_point(data = formatted_cov_dat, 
                       aes_string(x = "pos", y = "log2(cov + 1)", 
                                  alpha = fade_by, colour = color_by),
                       stroke = 0) +
            geom_smooth(data = formatted_cov_dat,
                        aes(x = pos, y = log2(cov + 1)),
                        se = FALSE, colour = myCbPal[2])
    } else {
        p_cov <- p_cov +
            geom_line(data = formatted_cov_dat,
                      aes_string(x = "pos", y = "log2(cov + 1)", 
                                 group = "lib_id", colour = color_by,
                                 alpha = fade_by),
                      stat = "smooth", method = "loess", 
                      se = FALSE, size = 1)
    }
    
    if (color_by == "lib_num") {
        p_cov <- p_cov + 
            scale_color_manual(values = cc) +
            guides(colour = FALSE)        
    } else {
        p_cov <- p_cov +
            scale_color_gradient(low = myCbPal[3], high = myCbPal[2])
    }
    
    if (is.numeric(fade_by)) {
        p_cov <- p_cov +
            scale_alpha_continuous(range = c(fade_by, fade_by)) +
            guides(alpha = FALSE)
    } else {
        p_cov <- p_cov +
            scale_alpha_continuous(range = c(0.2, 0.8))
    }
    
    if (color_by == fade_by) {
        p_cov <- p_cov +
            guides(alpha = FALSE)
    }
    
    p_cov <- p_cov +    
        scale_fill_viridis(discrete = TRUE) +
        theme_gray()
    
    if (split) {
        p_cov <- p_cov +
            facet_grid(donor_id ~ timepoint)
    }
    
    if (hide_legend) {
        p_cov <- p_cov +
            guides(fill = FALSE, colour = FALSE, alpha = FALSE)
    }
    return(p_cov)
}

# format data -------------------------------------------------------------
xcript_dat <- as.data.frame(xcripts_gtf) %>% 
    filter(seqnames != "NR_047551") # removing because non-coding RNA

# quick fix for egfrt_dat
egfrt_dat <- egfrt_dat %>% 
    dplyr::rename(egfr_xcript = transcript_id) %>% 
    filter(seqnames != "NR_047551") # removing because non-coding RNA


bulk_cov_dat <- bulk_lib_dat %>% 
    prep_cov_dat(car_cov_dat, bulk_metric_dat) %>% 
    filter(donor_id %in% c("x145", "x194", "x228"))

p89_c1_cov_dat <- sc_lib_dat %>% 
    prep_cov_dat(car_cov_dat, sc_metric_dat)

p85_c1_cov_dat <- p85_lib_dat %>% 
    prep_cov_dat(car_cov_dat, p85_metric_dat)



