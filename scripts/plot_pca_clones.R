
# load libs ---------------------------------------------------------------

library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggthemes)
library(viridis)

# load data ---------------------------------------------------------------

resids <- readRDS("data/from_masanao/zlm_resid_labeled.rds")
annotation <- readRDS("data/from_masanao/annotation.rds")

p89_c1_sample_dat <- read_csv("data/clean/p89_c1_compiled_sample_data.csv")

p89_c1_jxn_summary_dat <- read_csv("data/clean/p85_c1_tcr_summary_data.csv")
p89_c1_clonotype_dat <- read_csv("data/clean/p85_c1_compiled_tcr_data.csv")

cb_pal <- colorblind_pal()(8)

# load functions ----------------------------------------------------------

source("scripts/data_cleaning_functions.R")

# calculate principal components ------------------------------------------

pca <- prcomp(t(resids))

# format data -------------------------------------------------------------

pc_dat <- as.data.frame(pca$x) %>% 
    add_rownames(var = "wellKey") %>%
    select(wellKey, PC1, PC2) %>% 
    left_join(annotation %>% 
                  select(lib_id = libID, wellKey), .) %>% 
    select(-wellKey) %>% 
    inner_join(p89_c1_sample_dat %>% 
                   select(lib_id, donor_id, timepoint = timepoint_short) %>% 
                   left_join(p89_c1_clonotype_dat), .) %>% 
    replace_na(list(trav = "null", trbv = "null", clone_id = "null"))
    
# prep data for plotting --------------------------------------------------

# check for clones detected at multiple timepoints

# pc_dat <- pc_dat %>% 
#     group_by(donor_id, clone_id) %>% 
#     dplyr::mutate(num_timepoints = n_distinct(timepoint)) %>% 
#     ungroup() %>% 
#     mutate(jxn_detected = clone_id != "null",
#            num_timepoints = ifelse(jxn_detected, num_timepoints, 0)) %>% 
#     mutate(num_timepoints = factor(num_timepoints))

get_chain_status <- Vectorize(function(trav, trbv) {
    if(!(trav == "null") && !(trbv == "null")) {
        return("both")
    } else if((trav == "null") && (trbv == "null")) {
        return("none")
    } else if(!(trav == "null")) {
        return("alpha")
    } else {
        return("beta")
    }
})

# check for alleles detected at multiple timepoints
x <- pc_dat %>% 
    mutate(detected = get_chain_status(trav, trbv),
           any_detected = detected != "none") %>% 
    gather(chain, allele, trav:trbv) %>% 
    group_by(donor_id, chain, allele) %>% 
    mutate(num_timepoints = n_distinct(timepoint)) %>% 
    ungroup() %>% 
    mutate(num_timepoints = ifelse(allele != "null", num_timepoints, 0)) %>% 
    group_by(lib_id) %>% 
    mutate(repeat_alpha = ifelse(any(chain == "trav" & num_timepoints == 3), 
                                 allele[chain == "trav" & num_timepoints == 3], 
                                 "other"), 
           repeat_beta = ifelse(any(chain == "trbv" & num_timepoints == 3), 
                                allele[chain == "trbv" & num_timepoints == 3], 
                                "other"),
           num_timepoints = max(num_timepoints),
           any_repeated = any(repeat_alpha != "other") | 
               any(repeat_beta != "other"),
           highlight = any_detected & any_repeated) %>% 
    ungroup() %>% 
    spread(chain, allele)


# create plot -------------------------------------------------------------

# n_fill_colors <- x[["allele"]][x[["allele"]] != "null" & 
#                                    x[["chain"]] == "trbv"] %>% 
#     n_distinct()
n_fill_colors <- n_distinct(x[["trbv"]])
fill_cb_pal <- colorRampPalette(cb_pal)(n_fill_colors)
fill_cb_pal[1] <- "#444444"

# n_colour_colors <- x[["allele"]][x[["allele"]] != "null" & 
#                                      x[["chain"]] == "trav"] %>% 
#     n_distinct()
n_colour_colors <- n_distinct(x[["trav"]])
colour_cb_pal <- colorRampPalette(cb_pal)(n_colour_colors)
colour_cb_pal[1] <- "#444444"

# pc_dat %>% 
#     ggplot(aes(x = PC1, y = PC2)) +
#     geom_point(aes(fill = clone_id, alpha = jxn_detected, 
#                    size = num_timepoints, shape = num_timepoints)) +
#     facet_grid(donor_id ~ timepoint) +
#     scale_fill_manual(values = clone_cb_pal) +
#     scale_alpha_manual(values = c(0.2, 0.8)) +
#     scale_shape_manual(values = c(21, 24, 22, 23)) +
#     scale_size_manual(values = c(2, 2, 4, 6)) +
#     guides(fill = FALSE, size = FALSE) +
#     theme_bw() +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())

x %>% 
    ggplot(aes(x = PC1, y = PC2, 
               fill = trbv, colour = trav,
               alpha = any_detected, size = as.factor(num_timepoints))) +
    geom_point(stroke = 2) +
    geom_point(data = x %>% filter(detected == "none"), 
               shape = 21, stroke = 2) +
    geom_point(data = x %>% filter(detected == "alpha"), 
               shape = 24, stroke = 2) +
    geom_point(data = x %>% filter(detected == "beta"),
               shape = 22, stroke = 2) +
    geom_point(data = x %>% filter(detected == "both"),
               shape = 23, stroke = 2) +
    facet_grid(donor_id ~ timepoint) +
#     scale_fill_viridis(discrete = TRUE) +
#     scale_color_viridis(discrete = TRUE) +
    scale_fill_manual(values = fill_cb_pal) +
    scale_colour_manual(values = colour_cb_pal) +
    scale_alpha_manual(values = c(0.1, 0.7)) +
    scale_size_manual(values = c(1, 1, 2, 4)) +
    guides(colour = FALSE, #guide_legend(override.aes = list(shape = 21)),
           fill = FALSE, #guide_legend(override.aes = list(shape = 21, stroke = 0, size = 4)),
           size = guide_legend(override.aes = list(shape = 21)),
           alpha = guide_legend(override.aes = list(shape = 21))) +
    # scale_shape_manual(values = c(21, 24, 22, 23)) +
    theme_bw() # +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())

# other plot --------------------------------------------------------------

x %>% 
    gather(repeated_chain, repeated_allele, repeat_alpha:repeat_beta) %>% 
    filter(repeated_allele != "other") %>% 
    ggplot(aes(x = repeated_allele)) +
    geom_bar(aes(fill = donor_id), stat = "count") + 
    scale_fill_viridis(discrete = TRUE) + 
    facet_grid(~ timepoint)


# other plot 2 ------------------------------------------------------------

x %>% 
    filter(repeat_beta != "other") %>% 
    ggplot(aes(x = repeat_beta)) +
    geom_bar(aes(fill = donor_id), stat = "count") + 
    scale_fill_viridis(discrete = TRUE) + 
    facet_wrap(~ timepoint)

# get clones w/ multiple timepoints ---------------------------------------

clone_summary <- pc_dat %>% 
    filter(num_timepoints %in% c(2,3)) %>% 
    select(lib_id, donor_id, timepoint, clone_id, num_timepoints)

# save plot ---------------------------------------------------------------

# ggsave("car_detect_fig.png", combined_plot, 
#        scale = 1.6,
#        width = 17.35, height = 17.35, units = "cm", dpi = 300)
