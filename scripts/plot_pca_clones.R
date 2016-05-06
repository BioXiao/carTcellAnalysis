
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

# load("data/sample_metrics_data.RData")

p89_c1_sample_dat <- read_csv("data/clean/p89_c1_compiled_sample_data.csv")

# load("data/sample_tcr_data.RData")

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
pc_dat <- pc_dat %>% 
    group_by(donor_id, clone_id) %>% 
    dplyr::mutate(num_timepoints = n_distinct(timepoint)) %>% 
    ungroup() %>% 
    mutate(jxn_detected = clone_id != "null",
           num_timepoints = ifelse(jxn_detected, num_timepoints, 0)) %>% 
    mutate(num_timepoints = factor(num_timepoints))

# create plot -------------------------------------------------------------

n_colors <- n_distinct(pc_dat[["clone_id"]])
clone_cb_pal <- colorRampPalette(cb_pal)(n_colors)

pc_dat %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = clone_id, alpha = jxn_detected, 
                   size = num_timepoints, shape = num_timepoints)) +
    facet_grid(donor_id ~ timepoint) +
    scale_fill_manual(values = clone_cb_pal) +
    scale_alpha_manual(values = c(0.2, 0.8)) +
    scale_shape_manual(values = c(21, 24, 22, 23)) +
    scale_size_manual(values = c(2, 2, 4, 6)) +
    guides(fill = FALSE, size = FALSE) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

# get clones w/ multiple timepoints ---------------------------------------

clone_summary <- pc_dat %>% 
    filter(num_timepoints %in% c(2,3)) %>% 
    select(lib_id, donor_id, timepoint, clone_id, num_timepoints)

