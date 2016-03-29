
library(tximport)
library(readr)
library(dplyr)

xcript_list <- c("CAR-1", "NM_001561", "NM_000734", "NM_198053", "NM_001243077",
                 "NM_001243078", "NM_006139", "NM_000758", "NM_201283", 
                 "NM_201282", "NM_201284", "NM_005228", "NR_047551")

salmon_file_list <- list.files("~/Box Sync/data/projects/carTcellAnalysis/results/salmon/150910_D00565_0096_AC6REVANXX/P89-10",
    full.names = TRUE, recursive = TRUE) %>% 
    .[str_detect(., ".sf$")]

salmon_quant <- tximport(salmon_file_list[1], type = "salmon", txOut = TRUE,
                         reader = read_tsv)
salmon_counts <- salmon_quant$counts %>% 
    as.data.frame() %>% 
    add_rownames("tx_id") %>% 
    filter(tx_id %in% xcript_list)
