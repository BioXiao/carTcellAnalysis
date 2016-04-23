
# load libs ---------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(scales)
library(ggthemes)
library(viridis)

# load data ---------------------------------------------------------------
resids <- readRDS("data/zlm_resid_labeled.rds")
annotation <- readRDS("data/annotation.rds")

load("data/sample_metrics_data.RData")
load("data/sample_tcr_data.RData")

cb_pal <- colorblind_pal()(8)

# define functions --------------------------------------------------------

# clean up duplicated headers
clean_dup_names <- function(df) {
    df_names <- names(df)
    is_dup_name <- duplicated(df_names)
    df_names[is_dup_name] <- str_c(df_names[is_dup_name], "2")
    names(df) <- df_names
    return(df)
}

# simplify library ID
simplify_lib_id <- function(df) {
    df_names <- names(df)
    is_lib_name <- str_detect(tolower(df_names), "lib.*id")
    df_names[is_lib_name] <- "lib_id"
    names(df) <- df_names
    
    df %>% 
        mutate(lib_id = str_extract(lib_id, "lib[0-9]+"))
}

# clean/relabel donor ID
clean_donor_ids <- function(df) {
    df %>% 
        mutate(donor_id = tolower(donor_id))
}

# relabel timepoints
relabel_timepoints <- function(df) {
    df %>% 
        mutate(timepoint = str_replace(timepoint, " ", "")) %>% 
        mutate(timepoint = str_replace(timepoint, "InfusionProduct", "IP")) %>% 
        mutate(timepoint = str_replace(timepoint, "Day0", "t0")) %>% 
        mutate(timepoint = str_replace(timepoint, "Day(7|8|9|12)", "t1")) %>% 
        mutate(timepoint = str_replace(timepoint, "Day(26|28|29|33)", "t2"))
}

# format data -------------------------------------------------------------

# compute principal components
pca <- prcomp(t(resids))

pc_dat <- as.data.frame(pca$x) %>% 
    add_rownames(var = "wellKey") %>% 
    left_join(annotation %>% 
                  select(lib_id = libID, wellKey), .) %>% 
    select(-wellKey) %>% 
    inner_join(sc_lib_dat %>% 
                   clean_dup_names() %>% 
                   simplify_lib_id() %>%
                   clean_donor_ids() %>% 
                   select(lib_id, donor_id, timepoint) %>% 
                   relabel_timepoints() %>% 
                   left_join(jxn_dat), .) %>% 
    replace_na(list(TRAV = "null", TRBV = "null", clone_id = "null"))
    

# prep data for plotting --------------------------------------------------

# check for clones detected at multiple timepoints
pc_dat <- pc_dat %>% 
    group_by(donor_id, clone_id) %>% 
    dplyr::mutate(num_timepoints = n_distinct(timepoint)) %>% 
    ungroup() %>% 
    mutate(jxn_detected = clone_id != "null",
           num_timepoints = ifelse(jxn_detected, num_timepoints, 0)) %>% 
    mutate(num_timepoints = factor(num_timepoints))

# make plot ---------------------------------------------------------------
n_colors <- n_distinct(pc_dat$clone_id)
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

