
# load libs ---------------------------------------------------------------

library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(viridis)
library(cowplot)

# load data ---------------------------------------------------------------

p89_c1_sample_dat <- read_csv("data/clean/p89_c1_compiled_sample_data.csv")
p89_bulk_sample_dat <- read_csv("data/clean/p89_bulk_compiled_sample_data.csv")
p85_c1_sample_dat <- read_csv("data/clean/p85_c1_compiled_sample_data.csv")

car_cov_dat <- read_csv("data/clean/all_compiled_car_coverage_data.csv")
car_dat <- read_csv("data/clean/car_gtf_data.csv")

egfr_cov_dat <- read_csv("data/clean/all_compiled_egfr_coverage_data.csv")
egfrt_dat <- read_csv("data/clean/egfr_gtf_data.csv")

salmon_quant_dat <- read_csv("data/clean/all_compiled_abundance_data.csv")

lib_list <- readRDS("data/from_masanao/libs_NHL_after_filter.rds")

gff_file <- "data/annotation/_old/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)
xcript_dat <- as.data.frame(xcripts_gtf) %>% 
    mutate(seqnames = as.character(seqnames)) %>% 
    filter(seqnames != "NR_047551") # removing because non-coding RNA

# load functions ----------------------------------------------------------

source("scripts/data_cleaning_functions.R")

# define functions --------------------------------------------------------

# create/format data frame for plotting coverage
prep_cov_dat <- function(sample_dat, cov_dat, salmon_quant_dat, xcript_dat) {
    sample_dat %>% 
        select(lib_id, donor_id, timepoint) %>% 
        left_join(cov_dat) %>% 
        left_join(salmon_quant_dat %>% 
                      select(lib_id, name, tpm) %>% 
                      spread(name, tpm) %>% 
                      dplyr::rename(CAR = `CAR-1`) %>% 
                      relabel_transcripts(xcript_dat))
}

# crazy function to build plots
plot_coverage <- function(formatted_cov_dat, gtf_dat, 
                          fits_only = FALSE, all_fits = FALSE,
                          color_by = "lib_num") {
    
    v_pal <- viridis_pal()(8)
    
    # add lib numbers for coloring
    formatted_cov_dat <- formatted_cov_dat %>% 
        group_by(donor_id, timepoint) %>% 
        mutate(lib_num = dense_rank(lib_id) / n_distinct(lib_id),
               lib_num = as.character(lib_num)) %>% 
        ungroup()
    
    # determine plot height
    height <- log2(max(formatted_cov_dat[["cov"]], na.rm = TRUE) + 1)
    
    # build plot
    p_cov <- ggplot() +
        geom_rect(data = gtf_dat, 
                  aes(xmin = start, xmax = end, ymin = 0, ymax = height, 
                      fill = segment),
                  alpha = 0.5, colour = "white")
    if (!fits_only) {
        p_cov <- p_cov +
            geom_point(data = formatted_cov_dat, 
                       aes_string(x = "pos", y = "log2(cov + 1)"),
                       stroke = 0.5, size = 1.5, alpha = 0.3, 
                       colour = "#333333")
    }
    if (all_fits) {
        p_cov <- p_cov +
            geom_smooth(data = formatted_cov_dat,
                        aes_string(x = "pos", y = "log2(cov + 1)",
                                   group = "lib_id", colour = color_by),
                        method = "gam", formula = y ~ s(x), se = FALSE,
                        size = 1)
    } else {
        p_cov <- p_cov +
            geom_smooth(data = formatted_cov_dat,
                        aes(x = pos, y = log2(cov + 1)),
                        size = 1.5, colour = viridis_pal()(12)[11])
    }

    p_cov <- p_cov +
        scale_color_viridis(begin = 0, end = 1, option = "B") +
        scale_fill_viridis(discrete = TRUE, direction = -1) +
        theme_bw() +
        scale_x_continuous(expand = c(0.01, 0)) +
        scale_y_continuous(expand = c(0.01, 0), limits = c(0, height)) +
        theme(legend.position = "top")

    return(p_cov)
}

# even crazier function to assign CAR detection status based on coverage rules
detect_car <- function(formatted_cov_dat, cov_min = 1, run_min = 8) {
    cov_indicator_dat <- formatted_cov_dat %>% 
        select(lib_id, segment, pos, cov) %>% 
        group_by(lib_id) %>% 
        mutate(prev = lag(segment, 1), 
               subs = lead(segment, 1)) %>% 
        ungroup() %>% 
        mutate(junction = ifelse(!is.na(prev) & prev != subs, 
                                 str_c(prev, subs, sep = ":"), NA)) %>% 
        ungroup() %>% 
        filter(segment %in% c("CD19scFv", "T2A", "EGFRt") | !is.na(junction)) %>% 
        mutate(nz = cov >= cov_min) %>% 
        group_by(lib_id, segment) %>% 
        mutate(run_length = rep(rle(unlist(nz)) %>% .$lengths,
                                rle(unlist(nz)) %>% .$lengths)) %>% 
        ungroup() %>% 
        mutate(run_length = ifelse(nz, run_length, 0)) %>% 
        group_by(lib_id, junction) %>% 
        mutate(jxn_length = rep(rle(unlist(nz)) %>% .$lengths, 
                                rle(unlist(nz)) %>% .$lengths)) %>% 
        ungroup() %>% 
        mutate(jxn_length = ifelse(!is.na(junction) & nz, jxn_length, 0)) %>% 
        group_by(lib_id) %>% 
        summarise(car_detected = max(run_length) >= run_min |
                      sum(jxn_length == 2) >= 1)
    
    formatted_cov_dat %>% 
        left_join(cov_indicator_dat %>% 
                      select(lib_id, car_detected)) %>% 
        group_by(lib_id) %>% 
        mutate(nz_cov = sum(cov > 0) >= 5 | max(cov) >= 2) %>% 
        ungroup()
}

# format CAR coverage data -------------------------------------------------

p89_bulk_cov_dat <- p89_bulk_sample_dat %>% 
    prep_cov_dat(car_cov_dat, salmon_quant_dat, xcript_dat) %>% 
    filter(donor_id %in% c("x145", "x194", "x228"))

p89_c1_cov_dat <- p89_c1_sample_dat %>% 
    prep_cov_dat(car_cov_dat, salmon_quant_dat, xcript_dat) %>% 
    filter(lib_id %in% lib_list)

p85_c1_cov_dat <- p85_c1_sample_dat %>% 
    prep_cov_dat(car_cov_dat, salmon_quant_dat, xcript_dat)

# format EGFR coverage data -----------------------------------------------

p89_c1_egfr_cov_dat <- sc_lib_dat %>% 
    prep_cov_dat(egfr_cov_dat %>% 
                     filter(egfr_xcript != "EGFRt_r5"), # remove non-coding RNA
                 sc_metric_dat) %>% 
    filter(!is.na(egfr_xcript))

# assign CAR detection status in single cell libs -------------------------

p89_cov_w_car_detect_dat <- p89_c1_cov_dat %>% 
    detect_car(cov_min = 1, run_min = 8)

p89_cov_w_car_detect_dat %>% 
    filter(!is.na(nz_cov)) %>% 
    group_by(car_detected, nz_cov) %>% 
    summarise(n_libs = n_distinct(lib_id))

p85_cov_w_car_detect_dat <- p85_c1_cov_dat %>% 
    detect_car(cov_min = 1, run_min = 8)

p85_cov_w_car_detect_dat %>% 
    filter(!is.na(nz_cov)) %>% 
    group_by(car_detected, nz_cov) %>% 
    summarise(n_libs = n_distinct(lib_id))

# create CAR coverage plots -----------------------------------------------

car_plot <- bulk_cov_dat %>% 
    mutate(samples = "CAR T-cells (bulk)") %>% 
    bind_rows(p89_c1_cov_dat %>% 
                  mutate(samples = "CAR T-cells (single)")) %>% 
    bind_rows(p85_c1_cov_dat %>% 
                  mutate(samples = "MAI T-cells (single)")) %>% 
    plot_coverage(car_dat, color_by = "log2(CAR + 1)",
                  all_fits = FALSE, fits_only = FALSE) +
    guides(fill = guide_legend(title = "CAR segment", 
                               title.theme = element_text(size = 10, 
                                                          face = "bold", 
                                                          angle = 0),
                               nrow = 1)) +
    xlab("Position [bp]") +
    ylab("Coverage [log2(reads + 1)]") +
    facet_grid( ~ samples) +
    theme(panel.margin = unit(2, "lines"),
          plot.margin = unit(c(0, 1, 0, 1), "lines"))
    
# create EGFR coverage plots ----------------------------------------------

xcript_labels = list()
for (xcript in egfrt_dat$egfr_xcript) {
    xcript_labels[xcript] <- egfrt_dat$seqnames[which(egfrt_dat$egfr_xcript == xcript)]
}
xcript_labels <- unlist(xcript_labels)

egfr_plot <- p89_c1_egfr_cov_dat %>% 
    plot_coverage(egfrt_dat, all_fits = FALSE, fits_only = FALSE) +
    guides(fill = guide_legend(title = "CAR overlap", 
                               title.theme = element_text(size = 10, 
                                                          face = "bold", 
                                                          angle = 0),
                               nrow = 1)) +
    xlab("Position [bp]") +
    ylab("Coverage [log2(reads + 1)]") +
    facet_wrap(~ egfr_xcript, scales = "free_x",
               labeller = labeller(egfr_xcript = xcript_labels)) +
    theme(plot.margin = unit(c(0, 1, 0, 1), "lines"))

# create plot of CAR coverage for remaining libs ---------------------------

donor_id_labels <- c(x145 = "x145", x194 = "x194", x228 = "x228")
# donor_id_labels <- c(x145 = "ALL1", x194 = "ALL2", x228 = "ALL3")
timepoint_labels <- c(IP = "IP", t1 = "Peak", t2 = "Late")

car_plot_2 <- p89_cov_w_car_detect_dat %>% 
    filter(!is.na(nz_cov),
           car_detected) %>% 
    plot_coverage(car_dat, all_fits = TRUE, fits_only = TRUE,
                  color_by = "log2(CAR + 1)") +
    guides(fill = FALSE,
           colour = guide_colorbar(title = "CAR abundance [log2(TPM + 1)]", 
                                   title.theme = element_text(size = 10, 
                                                              face = "bold", 
                                                              angle = 0))) +
    xlab("Position [bp]") +
    ylab("Coverage [log2(reads + 1)]") +
    facet_grid(donor_id ~ timepoint, 
               labeller = labeller(donor_id = donor_id_labels,
                                   timepoint = timepoint_labels)) +
    theme(panel.margin = unit(1, "lines"),
          plot.margin = unit(c(0, 1, 0, 1), "lines"))

# build combined plot -----------------------------------------------------

combined_plot <- ggdraw() +
    draw_plot(car_plot, 0, 0.6, 1, 0.4) +
    draw_plot(egfr_plot, 0, 0, 0.45, 0.6) +
    draw_plot(car_plot_2, 0.45, 0, 0.55, 0.6) +
    draw_plot_label(c("A", "B", "C"), c(0, 0, 0.45), c(0.98, 0.58, 0.58), 
                    size = 12)

# save plot ---------------------------------------------------------------

ggsave("car_detect_fig.png", combined_plot, 
       scale = 1.6,
       width = 17.35, height = 17.35, units = "cm", dpi = 300)