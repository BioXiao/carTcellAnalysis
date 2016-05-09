
# load packages -----------------------------------------------------------

library(readr)
library(stringr)
library(tidyr)
library(dplyr)
library(purrr)

devtools::source_url("https://raw.githubusercontent.com/jaeddy/tcrtools/master/R/prep_junctions.R")

# define functions --------------------------------------------------------

# function to filter IMGT junctions
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

# load IMGT data ----------------------------------------------------------

imgt_file <- "data/from_hannah/NHL_IMGT_productive_TCRs_89_5_9_unique.txt"

# read and filter IMGT junctions
p89_c1_jxn_dat <- read_tsv(imgt_file) %>% 
    select(sample = libID,
           v_gene = `V.gene`,
           j_gene = `J.gene`,
           junction) %>% 
    filter_jxns() %>% 
    rename(lib_id = sample)

# create summary data -----------------------------------------------------

# summarize junction detection (old)
p89_c1_jxn_summary_dat <- p89_c1_jxn_dat %>%
    group_by(lib_id) %>% 
    summarise(tra_pos = sum(str_detect(v_gene, "TRA")),
           trb_pos = sum(str_detect(v_gene, "TRB"))) %>% 
    mutate(tcr_pos = ifelse(tra_pos + trb_pos >= 2, 1, 0))

# format clonotypes

# merge V-genes and junctions into a single string for each sample and chain
p89_c1_clonotype_dat <- p89_c1_jxn_dat %>% 
    mutate(jxn_chain = ifelse(str_detect(v_gene, "TRA"),
                          "alpha", "beta")) %>% 
    group_by(lib_id, jxn_chain) %>% 
    summarise(jxn_id = str_c(junction, collapse = "|"),
              v_id = str_c(v_gene, collapse = "|")) %>% 
    ungroup() %>% 
    complete(lib_id, jxn_chain, 
             fill = list(jxn_id = "null", v_id = "null"))

# combine alpha and beta chain ids to form clonotype ids
p89_c1_clonotype_dat <- p89_c1_clonotype_dat %>% 
    group_by(lib_id) %>% 
    mutate(clone_jxn_id = str_c(jxn_id, collapse = "_"),
           clone_v_id = str_c(v_id, collapse = "_")) %>% 
    ungroup() %>% 
    arrange(lib_id)


# save compiled files -----------------------------------------------------

write_csv(p89_c1_jxn_summary_dat, "data/clean/p85_c1_tcr_summary_data.csv")

write_csv(p89_c1_clonotype_dat, "data/clean/p85_c1_compiled_tcr_data.csv")

