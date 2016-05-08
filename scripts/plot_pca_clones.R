
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

# mark samples for which either an alpha or beta junction was detected; 'melt'
# alpha and beta junction columns into 'chain' and 'id' columns
x <- pc_dat %>% 
    rename(alpha = trav, beta = trbv) %>% 
    mutate(detected = get_chain_status(alpha, beta),
           any_detected = detected != "none") %>% 
    gather(jxn_chain, jxn_id, alpha:beta)

# count how many times (distinct timepoints) each unique alpha or beta junction 
# appears for each donor; set `n_timepoints` for null junctions to zero
x <- x %>% 
    group_by(donor_id, jxn_chain, jxn_id) %>% 
    mutate(num_timepoints = n_distinct(timepoint)) %>% 
    ungroup() %>% 
    mutate(num_timepoints = ifelse(jxn_id != "null", num_timepoints, 0)) 

# for each sample, store the junction id in a new column if it appears in all
# timepoints for a donor; because we'll just have one point per sample, keep
# the max number of timepoints across junctions (indicating the number of 
# time points that either a repeated alpha OR beta junction is detected for 
# that sample)
x <- x %>% 
    group_by(lib_id) %>% 
    mutate(repeat_alpha = ifelse(any(jxn_chain == "alpha" & num_timepoints == 3), 
                                 jxn_id[jxn_chain == "alpha" & 
                                            num_timepoints == 3], 
                                 "other"), 
           repeat_beta = ifelse(any(jxn_chain == "beta" & num_timepoints == 3), 
                                jxn_id[jxn_chain == "beta" & 
                                           num_timepoints == 3], 
                                "other"),
           num_timepoints = max(num_timepoints),
           any_repeated = any(repeat_alpha != "other") | 
               any(repeat_beta != "other")) %>% 
    ungroup()

# spread or 'cast' the previously melted chain/id info back into separate
# columns, so alpha and beta junction can be used independently for plotting
x <- x %>% 
    spread(jxn_chain, jxn_id) %>% 
    mutate(alpha = as.factor(alpha),
           beta = as.factor(beta))


# create plot -------------------------------------------------------------

# create fill palette with unique color for each beta junction (plus one 
# extra for non-detected junctions)
n_fill_colors <- n_distinct(x[["beta"]])
fill_cb_pal <- colorRampPalette(cb_pal)(n_fill_colors)
fill_cb_pal[1] <- "#FFFFFF"

# create color palette with unique color for each alpha junction (plus one 
# extra for non-detected junctions)
n_colour_colors <- n_distinct(x[["alpha"]])
colour_cb_pal <- colorRampPalette(cb_pal)(n_colour_colors)
colour_cb_pal[1] <- "#FFFFFF"

x %>% 
    # filter(any_repeated) %>% 
    ggplot(aes(x = PC1, y = PC2, 
               fill = beta, colour = alpha, label = clone_id,
               alpha = any_detected, size = as.factor(num_timepoints))) +
    geom_point(shape = 21, stroke = 2) +
    # geom_text() +
    facet_grid(donor_id ~ timepoint) +
    scale_fill_manual(values = fill_cb_pal) +
    scale_colour_manual(values = colour_cb_pal) +
    scale_alpha_manual(values = c(0.1, 0.7)) +
    scale_size_manual(values = c(1, 1, 2, 4)) +
    guides(colour = guide_legend(override.aes = list(size = 3)),
           fill = guide_legend(override.aes = list(stroke = 0, size = 4)),
           size = guide_legend(title = "times observed"),
           alpha = FALSE) +
    theme_gray() # +
#     theme(panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank())


# test --------------------------------------------------------------------

# get factor levels of repeated beta junctions
beta_repeated <- x %>% filter(num_timepoints == 3) %>% 
    select(beta) %>% 
    mutate(level_idx = as.numeric(beta)) %>% 
    distinct()

fill_cb_sub_pal <- fill_cb_pal[beta_repeated[["level_idx"]] %>% sort()]

# get factor levels of repeated beta junctions
alpha_repeated <- x %>% filter(num_timepoints == 3) %>% 
    select(alpha) %>% 
    mutate(level_idx = as.numeric(alpha)) %>% 
    distinct()

colour_cb_sub_pal <- colour_cb_pal[alpha_repeated[["level_idx"]] %>% sort()]

x %>% 
    filter(num_timepoints == 3) %>% 
    ggplot(aes(x = 1, y = clone_id, 
               fill = beta, colour = alpha, label = clone_id)) +
    geom_point(shape = 21, stroke = 2, size = 4) +
    facet_wrap(~ donor_id, nrow = 3, scales = "free") +
    scale_fill_manual(values = fill_cb_sub_pal) +
    scale_colour_manual(values = colour_cb_sub_pal) +
    guides(colour = FALSE,
           fill = FALSE) +
    theme_gray()

# other plot --------------------------------------------------------------

x %>% 
    gather(repeated_chain, repeated_allele, repeat_alpha:repeat_beta) %>% 
    filter(repeated_allele != "other") %>% 
    ggplot(aes(x = repeated_allele)) +
    geom_bar(aes(fill = timepoint), stat = "count") + 
    scale_fill_viridis(discrete = TRUE) + 
    # coord_flip() +
    facet_wrap(~ donor_id, nrow = 3) +
    coord_flip()


# other plot 2 ------------------------------------------------------------

x %>% 
    filter(repeat_beta != "other") %>% 
    ggplot(aes(x = repeat_beta)) +
    geom_bar(aes(fill = donor_id), stat = "count") + 
    scale_fill_viridis(discrete = TRUE) + 
    facet_wrap(~ donor_id)

# get clones w/ multiple timepoints ---------------------------------------

clone_summary <- pc_dat %>% 
    filter(num_timepoints %in% c(2,3)) %>% 
    select(lib_id, donor_id, timepoint, clone_id, num_timepoints)

# save plot ---------------------------------------------------------------

# ggsave("car_detect_fig.png", combined_plot, 
#        scale = 1.6,
#        width = 17.35, height = 17.35, units = "cm", dpi = 300)
