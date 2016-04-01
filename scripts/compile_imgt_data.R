
library(dplyr)
library(stringr)
library(tidyr)
library(devtools)
jxn_script_url <- "https://raw.githubusercontent.com/jaeddy/tcrSeqAnalysis/master/R/prep_junctions.R"
source_url(jxn_script_url)

imgt_file <- "data/NHL_IMGT_productive_TCRs_89_5_9_unique.txt"
# Read and filter IMGT junctions
jxns <- read.delim(imgt_file, stringsAsFactors = FALSE) %>% 
    select(lib_id = libID,
           v_gene = `V.gene`,
           j_gene = `J.gene`,
           junction) %>% 
    filter_imgt_jxns()

# Summarize junction detection
jxn_summary_dat <- jxns %>% 
    group_by(lib_id) %>% 
    summarise(tra_pos = sum(str_detect(v_gene, "TRA")),
           trb_pos = sum(str_detect(v_gene, "TRB"))) %>% 
    mutate(tcr_pos = ifelse(tra_pos + trb_pos >= 2, 1, 0))

# Save image
save(jxn_summary_dat, file = "data/sample_tcr_data.RData")