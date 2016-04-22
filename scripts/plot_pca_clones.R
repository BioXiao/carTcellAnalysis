lib_list <- readRDS("data/libs_NHL_after_filter.rds")
resids <- readRDS("data/zlm_resid.rds")

load("data/sample_metrics_data.RData")
load("data/sample_tcr_data.RData")


library(dplyr)
library(ggplot2)
library(stringr)
library(scales)
library(ggthemes)
library(viridis)

myCbPal <- colorblind_pal()(8)
myCbPal[c(1, 3, 6)] <- myCbPal[c(6, 1, 3)]
myCbPal[3] <- "#666666"

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

# Compute principal components
pca <- prcomp(t(resids))

# Plot first two PCs and color by group
ggbiplot(pca, choices = c(1,2), var.axes = FALSE) + , groups = groupDat$group) +
    scale_color_colorblind() +
    guides(size = FALSE) +
    theme_classic() +
    theme(axis.title = element_text(size = 11))

pc_dat <- as.data.frame(pca$x) %>% 
    mutate(lib_id = lib_list) %>% 
    inner_join(sc_lib_dat %>% 
                   clean_dup_names() %>% 
                   simplify_lib_id() %>%
                   clean_donor_ids() %>% 
                   select(lib_id, donor_id, timepoint) %>% 
                   relabel_timepoints(), .)

head(pc_dat)[, 1:5] %>% 
    mutate(timepoint = str_replace(timepoint, " ", "")) %>% 
    mutate(timepoint = str_replace(timepoint, "InfusionProduct", "IP"))

pc_dat %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point() +
    geom_density_2d() +
    facet_grid(donor_id ~ timepoint)


construct_tcrs <- function(jxns_df, any = FALSE) {
    
    trav_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRAV")) %>% 
        select(-j_gene) %>% 
        dplyr::rename(trav_gene = v_gene, 
                      trav_jxn = junction)
    
    trbv_df <- jxns_df %>% 
        filter(str_detect(v_gene, "TRBV")) %>% 
        select(-j_gene) %>% 
        dplyr::rename(trbv_gene = v_gene,
                      trbv_jxn = junction)
    
    if (!any) {
        tcr_df <- inner_join(trav_df, trbv_df, by = "lib_id")
    } else {
        tcr_df <- full_join(trav_df, trbv_df, by = "lib_id") %>% 
            mutate(trav_gene = ifelse(is.na(trav_gene), "no_trav_gene", trav_gene),
                   trbv_gene = ifelse(is.na(trbv_gene), "no_trbv_gene", trbv_gene),
                   trav_jxn = ifelse(is.na(trav_jxn), "no_trav_jxn", trav_jxn),
                   trbv_jxn = ifelse(is.na(trbv_jxn), "no_trbv_jxn", trbv_jxn))
    }
    return(tcr_df)
}

tcr_dat <- jxns %>% 
    construct_tcrs(any = TRUE) %>% 
    mutate(tcr_str = str_c(trav_gene, trbv_gene, sep = ":"))

pc_dat <- pc_dat %>% 
    left_join(tcr_dat %>% 
                  select(lib_id, trav_gene, trbv_gene, tcr_str))




# make plot ---------------------------------------------------------------

color_by <- "trav_gene"
fill_by <- "trbv_gene"
num_colors <- n_distinct(pc_dat[[color_by]]) - 1
# cc <- colorRampPalette(c(myCbPal[1], myCbPal[2], myCbPal[4]))(num_colors)
color_pal <- c("#000000", colorRampPalette(myCbPal)(num_colors))
# color_pal <- viridis_pal()(num_colors)

num_fills <- n_distinct(pc_dat[[fill_by]]) - 1
fill_pal <- c("#000000", colorRampPalette(myCbPal)(num_fills))
# fill_pal <- viridis_pal()(num_colors)

pc_dat %>% 
#     filter(trav_gene != "no_trav_gene",
#            trbv_gene != "no_trbv_gene") %>% 
    ggplot(aes(x = PC1, y = PC2)) +
    geom_point(aes_string(colour = color_by, fill = fill_by), 
               shape = 21, size = 4, alpha = 0.8) +
    # geom_density_2d() +
    scale_color_manual(values = color_pal) +
    scale_fill_manual(values = fill_pal) +
    facet_grid(donor_id ~ timepoint)
