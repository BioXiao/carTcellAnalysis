
library(readr)
library(stringr)
library(dplyr)

xcript_list <- c("CAR-1", "NM_001561", "NM_000734", "NM_198053", "NM_001243077",
                 "NM_001243078", "NM_006139", "NM_000758", "NM_201283", 
                 "NM_201282", "NM_201284", "NM_005228", "NR_047551")

build_quant_dat <- function(salmon_files) {
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
    lapply(as.list(salmon_files), function(x) {
        lib_id <- str_extract(x, "lib[0-9]+")
        map_rate <- read_lines(x, skip = 7, n_max = 1) %>% 
            str_extract("[0-9]+\\.[0-9]+") %>% 
            as.numeric() / 100
        data_frame(lib_id = lib_id, map_rate = map_rate)
    }) %>% 
        bind_rows()
}

# collect Salmon data for all samples
salmon_files <- list.files("~/Box Sync/data/projects/carTcellAnalysis/results/salmon",
                           full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".sf$")]

salmon_quant_dat <- build_quant_dat(salmon_files)
map_rate_dat <- build_map_rate_dat(salmon_files)

# save image
save(salmon_quant_dat, map_rate_dat, file = "data/sample_salmon_data.RData")
