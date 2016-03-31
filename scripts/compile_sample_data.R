
# load packages -----------------------------------------------------------

library(readr)
library(readxl)
library(stringr)
library(dplyr)
library(tidyr)


# specify file paths ------------------------------------------------------

sample_file_list <- list(
    p89_5_ip = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-5/P89-5_LibAnnotation_x194IP_2_1771023248_24Aug2015.xlsx",
    p89_5_d7 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-5/P89-5_C1 LibAnnotaion_x194Day7_1771023246_19Aug2015.xlsx",
    p89_5_d28 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-5/P89-5_C1 LibAnnotation_x194Day28_1771029069_2015Jul27.xlsx",
    p89_7_ip = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-7/P89-7_C1 LibAnnotation_x145_IP_1771023244_2015Aug12.xlsx",
    p89_7_d8 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-7/P89-7_2015Jul16_C1 LibAnnotation_1771023223_x145Day8_1.xlsx",
    p89_7_d26 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-7/P89-7_C1 LibAnnotation_x145Day26_1771023239_04Aug2015.xlsx",
    p89_9_ip = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-9/P89-9_C1 LibAnnotation_1771023229_x228IP_2015Jul23.xlsx",
    p89_9_d7 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-9/P89-9_2015Jul21_C1 LibAnnotation_x228Day7_1771023227.xlsx",
    p89_9_d28 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/P89-9/P89-9_C1 LibAnnotation_x228_Day28_1771023245_2015Aug18.xlsx"
) 

p89_bulk <- "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/For analyst/P89 - Cameron/Bulk samples P89-6, -8, -10. -11, -12, -13/P89 bulk libs.xlsx"

p85_sample_file_list <- list(
    p85_5_1 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-5 C1 LibAnnotation (1of6) 1771023212_12May2015.xlsx",
    p85_5_2 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-5 C1 LibAnnotation (2of6) 1771023213_13May2015.xlsx",
    p85_5_3 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-5 C1 LibAnnotation (3of6) 1771023214_14May2015.xlsx",
    p85_6_1 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-6 C1 LibAnnotation (4of6) 1771023215_18May2015.xlsx",
    p85_6_2 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-6 C1 LibAnnotation (5of6) 1771023216_19May2015.xlsx",
    p85_6_3 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/For Analyst/P85 Prlic/P85-6 C1 LibAnnotation (6of6) 1771023217_20May2015.xlsx"
)

# specify metric file paths -----------------------------------------------

metric_file_list <- list(
    p89_5 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-5Processed_150916/metrics/P89-5_C6REVANXX_150916_combined_metrics.csv",
    p89_7_1 = "/Volumes/genomics/Illumina/150729_D00565_0092_AC6VC6ANXX/Project_P89-7Processed_150916/metrics/P89-7_C6VC6ANXX_150916_combined_metrics.csv",
    p89_7_2 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-7Processed_150916/metrics/P89-7_C6REVANXX_150916_combined_metrics.csv",
    p89_9_1 = "/Volumes/genomics/Illumina/150729_D00565_0092_AC6VC6ANXX/Project_P89-9Processed_150916/metrics/P89-9_C6VC6ANXX_150916_combined_metrics.csv",
    p89_9_2 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-9Processed_150916/metrics/P89-9_C6REVANXX_150916_combined_metrics.csv"
)

bulk_metric_file_list <- list(
    p89_6 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-6Processed_150917/metrics/P89-6_C6REVANXX_150917_combined_metrics.csv",
    p89_8 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-8Processed_150917/metrics/P89-8_C6REVANXX_150917_combined_metrics.csv",
    p89_10 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-10Processed_150917/metrics/P89-10_C6REVANXX_150917_combined_metrics.csv",
    p89_11 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-11Processed_150917/metrics/P89-11_C6REVANXX_150917_combined_metrics.csv",
    p89_12 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-12Processed_150917/metrics/P89-12_C6REVANXX_150917_combined_metrics.csv",
    p89_13 = "/Volumes/genomics/Illumina/150910_D00565_0096_AC6REVANXX/Project_P89-13Processed_150917/metrics/P89-13_C6REVANXX_150917_combined_metrics.csv"
)

p85_metric_file_list <- list(
    p85_5 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/Project_P85-5Processed_new/metrics/C6VC8ANXX_combined_metrics.csv",
    p85_6 = "/Volumes/genomics/Illumina/150528_D00565_0086_BC6VC8ANXX/Project_P85-6Processed_new/metrics/C6VC8ANXX_combined_metrics.csv"
)
# utility functions -------------------------------------------------------

snake_case <- function(str) {
    str_to_lower(str) %>% 
        str_replace_all(' ', '_')
}

df_to_snake_case <- function(df) {
    names(df) <- snake_case(names(df))
    return(df)
}
# load sample data ---------------------------------------------------------------

sc_lib_dat <- lapply(sample_file_list, function(x) {
    read_excel(x) %>% 
        df_to_snake_case()
}) %>% 
    Reduce(full_join, .) %>% 
    filter(!is.na(libid))

bulk_lib_dat <- read_excel(p89_bulk) %>% 
    .[, -ncol(.)] %>% 
    df_to_snake_case() %>% 
    separate(sample_name, c("donor_id", "timepoint"), 
             remove = FALSE, sep = "_", extra = "drop")

p85_lib_dat <- lapply(p85_sample_file_list, function(x) {
    read_excel(x) %>% 
        df_to_snake_case()
}) %>% 
    Reduce(full_join, .) %>% 
    filter(!is.na(libid))

# load metrics data -------------------------------------------------------

sc_metric_dat <- lapply(metric_file_list, function(x) {
    read_csv(x) %>% 
        df_to_snake_case()
}) %>% 
    bind_rows() %>% 
    filter(!is.na(libid))

bulk_metric_dat <- lapply(bulk_metric_file_list, function(x) {
    read_csv(x) %>% 
        df_to_snake_case()
}) %>% 
    bind_rows() %>% 
    filter(!is.na(libid))

p85_metric_dat <- lapply(p85_metric_file_list, function(x) {
    read_csv(x, na = "?") %>% 
        df_to_snake_case()
}) %>% 
    bind_rows() %>% 
    filter(!is.na(libid))

# save image
save(bulk_lib_dat, bulk_metric_dat, 
     sc_lib_dat, sc_metric_dat, 
     p85_lib_dat, p85_metric_dat,
     file = "data/sample_metrics_data.RData")
