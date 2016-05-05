
# load packages -----------------------------------------------------------

library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)
library(parallel)

# load functions ----------------------------------------------------------

source("scripts/data_cleaning_functions.R")

# specify file paths ------------------------------------------------------

p89_c1_sample_file_list <- list.files("data/sample_data", pattern = "P89", 
                                      full.names = TRUE) %>% 
    .[!str_detect(., "bulk")] %>% 
    as.list()

p89_bulk_sample_file_list <- list.files("data/sample_data", pattern = "P89", 
                                        full.names = TRUE) %>% 
    .[str_detect(., "bulk")] %>% 
    as.list()

p85_c1_sample_file_list <- list.files("data/sample_data", pattern = "P85",
                                      full.names = TRUE) %>% 
    as.list()

# specify metric file paths -----------------------------------------------

p89_c1_metric_file_list <- list.files("data/seq_metrics", pattern = "P89", 
                                      full.names = TRUE) %>% 
    .[str_detect(., "P89-(5|7|9)")] %>% 
    as.list()

p89_bulk_metric_file_list <- list.files("data/seq_metrics", pattern = "P89", 
                                        full.names = TRUE) %>% 
    .[!str_detect(., "P89-(5|7|9)")] %>% 
    as.list()

p85_c1_metric_file_list <- list.files("data/seq_metrics", pattern = "P85", 
                                      full.names = TRUE) %>% 
    as.list()

# load sample data ---------------------------------------------------------------

load_sample_data <- function(file_list) {
    sample_dat <- mclapply(file_list, function(x) {
        df <- read_excel(x) %>% 
            clean_headers() %>% 
            clean_dup_names() %>%
            select(matches(".+")) %>% 
            simplify_lib_id() 
        
        if(!("donor_id" %in% names(df))) {
            df <- df %>% 
                separate(sample_name, c("donor_id", "timepoint"), 
                         remove = FALSE, sep = "_", extra = "drop")
        }
        
        df %>% 
            relabel_timepoints() %>% 
            clean_donor_ids()
    }) %>% 
        bind_rows() %>% 
        filter(!is.na(lib_id))
    return(sample_dat)
}

p89_c1_sample_dat <- load_sample_data(p89_c1_sample_file_list)

p89_bulk_sample_dat <- load_sample_data(p89_bulk_sample_file_list)

p85_c1_sample_dat <- load_sample_data(p85_c1_sample_file_list)

# load metrics data -------------------------------------------------------

load_metric_data <- function(file_list) {
    metric_dat <- mclapply(file_list, function(x) {
        read_csv(x, na = c("", "?")) %>% 
            clean_headers() %>% 
            clean_dup_names() %>% 
            simplify_lib_id()
    }) %>% 
        bind_rows() %>% 
        filter(!is.na(lib_id))
}

p89_c1_metric_dat <- load_metric_data(p89_c1_metric_file_list)

p89_bulk_metric_dat <- load_metric_data(p89_bulk_metric_file_list)

p85_c1_metric_dat <- load_metric_data(p85_c1_metric_file_list)

# save compiled files -----------------------------------------------------

write_csv(p89_c1_sample_dat, "data/sample_data/p89_c1_compiled_sample_data.csv")
write_csv(p89_bulk_sample_dat, "data/sample_data/p89_bulk_compiled_sample_data.csv")
write_csv(p85_c1_sample_dat, "data/sample_data/p89_c1_compiled_sample_data.csv")

write_csv(p89_c1_metric_dat, "data/seq_metrics/p89_c1_compiled_metric_data.csv")
write_csv(p89_bulk_metric_dat, "data/seq_metrics/p89_bulk_compiled_metric_data.csv")
write_csv(p85_c1_metric_dat, "data/seq_metrics/p89_c1_compiled_metric_data.csv")

