
library(dplyr)
library(stringr)
library(tidyr)

# Function to filter IMGT junctions
filter_imgt_jxns <- function(imgt_jxns, min_length = 6) {
    
    imgt_jxns <- imgt_jxns %>% 
        filter(str_detect(v_gene, "^((?![C-G]).)*$"),
               str_detect(j_gene, "^((?![C-G]).)*$"),
               str_detect(junction, "^C"),
               str_detect(junction, "^((?!(\\*|_)).)*$"),
               !duplicated(.[, c(1:4)]),
               str_length(junction) > min_length)
    
    return(imgt_jxns)
}

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