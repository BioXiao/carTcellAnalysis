
# Load packages -----------------------------------------------------------

library(ggplot2)
library(readr)
library(stringr)
library(reshape2)
library(dplyr)


# Locate files ------------------------------------------------------------

# assumes result files located in the following subfolder
result_folder <- "data/results"

# # get list of all low confidence read count files for relevant flowcells
# gene_hits_files <- list.files(result_folder, full.names = TRUE) %>% 
#     .[str_detect(., "150(729|910).csv")]
# 
# # get list of all high confidence read count files for relevant flowcells
# gene_hits_hc_files <- list.files(result_folder, full.names = TRUE) %>% 
#     .[str_detect(., "150(729|910)_high_conf.csv")]

# for P85
# get list of all low confidence read count files for relevant flowcells
gene_hits_files <- list.files(result_folder, full.names = TRUE) %>% 
    .[str_detect(., "150528.csv")]

# get list of all high confidence read count files for relevant flowcells
gene_hits_hc_files <- list.files(result_folder, full.names = TRUE) %>% 
    .[str_detect(., "150528_high_conf.csv")]

# Load data ---------------------------------------------------------------

# read and merge all low confidence read count files
gene_hits_dat <- lapply(as.list(gene_hits_files), function(x) {
    df <- read_csv(x)
    names(df)[1] <- "lib_id"
    return(df)
}) %>% 
    bind_rows()

# read and merge all high confidence read count files
gene_hits_hc_dat <- lapply(as.list(gene_hits_hc_files), function(x) {
    df <- read_csv(x)
    names(df)[1] <- "lib_id"
    return(df)
}) %>% 
    bind_rows()


# Join, score data --------------------------------------------------------

# combine low and high confidence read counts, score libs as weak, medium, or
# strong supporting evidence (thresholds are somewhat arbitrary):
#   ~ weak - count > 0 in low confidence data
#   ~ medium - count > median count in low confidence data
#   ~ strong - count > 0 in high confidence data
gene_hits_all <- gene_hits_dat %>%
    select(-tpm) %>% 
    left_join(gene_hits_hc_dat %>% 
                  select(-tpm) %>% 
                  rename(num_reads_hc = num_reads)) %>% 
    mutate(weak = ifelse(num_reads > 0, TRUE, FALSE),
           medium = ifelse(num_reads > median(num_reads), TRUE, FALSE),
           strong = ifelse(num_reads_hc > 0, TRUE, FALSE))


# Summarize data ----------------------------------------------------------

# overall rates of CAR gene detection with weak, medium, or strong evidence
gene_hits_all %>% 
    summarise_each_(funs(mean), vars = c("weak", "medium", "strong"))

# rates by project
gene_hits_all %>% 
    group_by(project) %>% 
    summarise_each_(funs(mean), vars = c("weak", "medium", "strong"))

# rates by project and flowcell
gene_hits_all %>% 
    group_by(flowcell, project) %>% 
    summarise_each_(funs(mean), vars = c("weak", "medium", "strong"))


# Plot gene counts --------------------------------------------------------

# unwieldy bar graph of CAR gene read counts per sample
gene_hits_all %>%
    filter(medium) %>% 
    melt(measure.vars = c("num_reads", "num_reads_hc"),
         variable.name = "confidence", value.name = "num_reads") %>% 
    mutate(confidence = ifelse(str_detect(confidence, "hc"), 
                               "high", "low")) %>% 
    ggplot(aes(x = lib_id, y = num_reads)) +
    geom_bar(aes(fill = confidence), 
             stat = "identity", position = "dodge") +
    coord_flip()

