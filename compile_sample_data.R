
# load packages -----------------------------------------------------------

library(readr)
library(readxl)
library(stringr)
library(dplyr)


# specify file paths ------------------------------------------------------

file_list <- list(
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


# utility functions -------------------------------------------------------

snake_case <- function(str) {
    str_to_lower(str) %>% 
        str_replace_all(' ', '_')
}

df_to_snake_case <- function(df) {
    names(df) <- snake_case(names(df))
    return(df)
}
# load data ---------------------------------------------------------------

sc_lib_dat <- lapply(file_list, function(x) {
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

