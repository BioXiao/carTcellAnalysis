
# load libs ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(rtracklayer)
library(ggthemes)
library(viridis)
library(cowplot)


# load data ---------------------------------------------------------------

load("data/sample_metrics_data.RData")
load("data/sample_rapmap_data.RData")
load("data/sample_salmon_data.RData")
lib_list <- readRDS("data/libs_NHL_after_filter.rds")

gff_file <- "data/annotation/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)

source("scripts/data_cleaning_functions.R")

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
                          split = TRUE, fits_only = FALSE, all_fits = FALSE,
                          color_by = "lib_num", fade_by = "log2(CAR + 1)",
                          hide_legend = FALSE) {
    
    cb_pal <- colorblind_pal()(8)
    v_pal <- viridis_pal()(8)
    
    # add lib numbers for coloring
    formatted_cov_dat <- formatted_cov_dat %>% 
        group_by(donor_id, timepoint) %>% 
        mutate(lib_num = dense_rank(lib_id) / n_distinct(lib_id),
               lib_num = as.character(lib_num)) %>% 
        ungroup()
    
    # add even/odd indicators for label positioning
#     gtf_dat <- gtf_dat %>% 
#         mutate(seg_adjust = ifelse(row_number() %% 2, 0, -0.2),
#                seg_adjust = ifelse(row_number() %% 3, seg_adjust, -0.4))
#     print(gtf_dat)
#     
    # determine plot height
    height <- log2(max(formatted_cov_dat$cov, na.rm = TRUE) + 1)
    
    # build plot
    p_cov <- ggplot() +
        geom_rect(data = gtf_dat, 
                  aes(xmin = start, xmax = end, ymin = 0, ymax = height, 
                      fill = segment),
                  alpha = 0.5, colour = "white")
#         geom_label(data = gtf_dat,
#                   aes(x = start, y = height + seg_adjust, label = segment),
#                   hjust = "inward", vjust = "inward", size = 3,
#                   fill = "white")
    if (!fits_only) {
        p_cov <- p_cov +
            geom_point(data = formatted_cov_dat, 
                       aes_string(x = "pos", y = "log2(cov + 1)"),
                       stroke = 0.5, size = 2.5, alpha = 0.7, 
                       colour = "slategray")
    }
    if (all_fits) {
        p_cov <- p_cov +
#             geom_line(data = formatted_cov_dat,
#                       aes_string(x = "pos", y = "log2(cov + 1)", 
#                                  group = "lib_id", colour = color_by),
#                       stat = "smooth", method = "gam", formula = y ~ s(x),
#                       size = 1.5)
            geom_smooth(data = formatted_cov_dat,
                        aes_string(x = "pos", y = "log2(cov + 1)",
                                   group = "lib_id", colour = color_by),
                        method = "gam", formula = y ~ s(x), se = FALSE,
                        size = 1.5)
    } else {
        p_cov <- p_cov +
#             geom_line(data = formatted_cov_dat,
#                       aes(x = pos, y = log2(cov + 1)),
#                       stat = "smooth", method = "gam", formula = y ~ s(x), span = 0.5,
#                       size = 4, alpha = 0.5,
#                       colour = "slategray") +
#             geom_line(data = formatted_cov_dat,
#                       aes(x = pos, y = log2(cov + 1)),
#                       stat = "smooth", method = "gam", formula = y ~ s(x),
#                       size = 2, alpha = 1,
#                       colour = viridis_pal()(12)[11]) +
            geom_smooth(data = formatted_cov_dat,
                        aes(x = pos, y = log2(cov + 1)),
                        size = 2, colour = viridis_pal()(12)[11])
    }


    if (color_by == "lib_num") {
        p_cov <- p_cov + 
            scale_color_manual(values = line_cols) +
            guides(colour = FALSE)        
    } else {
        p_cov <- p_cov +
            scale_color_viridis(begin = 0, end = 1, option = "B")
            # scale_color_gradient(low = v_pal[1], high = v_pal[3])
    }
    
    p_cov <- p_cov +    
        scale_fill_viridis(discrete = TRUE, direction = -1) +
        theme_bw() +
        scale_x_continuous(expand = c(0.01, 0)) +
        scale_y_continuous(expand = c(0.01, 0), limits = c(0, height)) +
        theme(legend.position = "top")

    if (split) {
        p_cov <- p_cov +
            facet_grid(donor_id ~ timepoint)
    }
    
#     if (hide_legend) {
#         p_cov <- p_cov +
#             guides(fill = FALSE, colour = FALSE, alpha = FALSE)
#     }
    return(p_cov)
}

plot_coverage(bulk_cov_dat, car_dat, color_by = "log2(CAR + 1)", 
              all_fits = FALSE, fits_only = FALSE, split = FALSE)

# format reference data ----------------------------------------------------
xcript_dat <- as.data.frame(xcripts_gtf) %>% 
    filter(seqnames != "NR_047551") # removing because non-coding RNA

# quick fix for egfrt_dat
egfrt_dat <- egfrt_dat %>% 
    dplyr::rename(egfr_xcript = transcript_id) %>% 
    filter(seqnames != "NR_047551") # removing because non-coding RNA


# format coverage data ----------------------------------------------------

bulk_cov_dat <- bulk_lib_dat %>% 
    prep_cov_dat(car_cov_dat, bulk_metric_dat) %>% 
    filter(donor_id %in% c("x145", "x194", "x228"))

p89_c1_cov_dat <- sc_lib_dat %>% 
    prep_cov_dat(car_cov_dat, sc_metric_dat)

p85_c1_cov_dat <- p85_lib_dat %>% 
    prep_cov_dat(car_cov_dat, p85_metric_dat)


# create combined plots ---------------------------------------------------

bulk_cov_dat %>% 
    mutate(samples = "CAR T-cells (bulk)") %>% 
    bind_rows(p89_c1_cov_dat %>% 
                  mutate(samples = "CAR T-cells (single)")) %>% 
    bind_rows(p85_c1_cov_dat %>% 
                  mutate(samples = "MAI T-cells (single)")) %>% 
    plot_coverage(car_dat, color_by = "log2(CAR + 1)",
                  all_fits = FALSE, fits_only = FALSE, split = FALSE) +
    facet_grid( ~ samples)

# plot_grid(p1, p2, p3, ncol = 3, nrow = 1, 
#           hjust = 0)
