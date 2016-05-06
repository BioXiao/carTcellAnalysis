
# load packages -----------------------------------------------------------

library(rtracklayer)
library(readr)
library(stringr)
library(dplyr)
library(parallel)

# load data ---------------------------------------------------------------

gff_file <- "data/annotation/_old/carPlusRef.gtf"
xcripts_gtf <- import.gff2(gff_file)

# load functions ----------------------------------------------------------

source("scripts/data_cleaning_functions.R")

# define functions --------------------------------------------------------

build_quant_dat <- function(salmon_files, xcript_list) {
    lapply(as.list(salmon_files), function(x) {
        lib_id <- str_extract(x, "lib[0-9]+")
        salmon_quant <- read_tsv(x, skip = 8) %>% 
            rename(Name = `# Name`) %>% 
            filter(Name %in% xcript_list) %>% 
            mutate(lib_id = lib_id)
    }) %>% 
        bind_rows()
}

build_map_rate_dat <- function(salmon_files) {
    mclapply(as.list(salmon_files), function(x) {
        lib_id <- str_extract(x, "lib[0-9]+")
        map_rate <- read_lines(x, skip = 7, n_max = 1) %>% 
            str_extract("[0-9]+\\.[0-9]+") %>% 
            as.numeric() / 100
        data_frame(lib_id = lib_id, map_rate = map_rate)
    }) %>% 
        bind_rows()
}

# get transcript IDs ------------------------------------------------------

# extract IDs of transcripts in 'xcripts_gtf'
xcript_list <- xcripts_gtf %>% 
    seqnames() %>% 
    as.list() %>% 
    unique()

# specify file paths ------------------------------------------------------

salmon_files <- list.files("data/_salmon_old",
                           full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".sf$")]

# collect Salmon abundance data for all samples ---------------------------

salmon_quant_dat <- build_quant_dat(salmon_files, xcript_list) %>% 
    clean_headers()

# collect Salmon map rate data for all samples ----------------------------

map_rate_dat <- build_map_rate_dat(salmon_files)

# save compiled files -----------------------------------------------------

write_csv(salmon_quant_dat, "data/clean/all_compiled_abundance_data.csv")
write_csv(map_rate_dat, "data/clean/all_compiled_map_rate_data.csv")
