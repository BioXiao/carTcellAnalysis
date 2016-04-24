
# clean up duplicated headers
clean_dup_names <- function(df) {
    df_names <- names(df)
    is_dup_name <- duplicated(df_names)
    df_names[is_dup_name] <- str_c(df_names[is_dup_name], "2")
    names(df) <- df_names
    return(df)
}

# simplify library ID
simplify_lib_id <- function(df) {
    df_names <- names(df)
    is_lib_name <- str_detect(tolower(df_names), "lib.*id")
    df_names[is_lib_name] <- "lib_id"
    names(df) <- df_names
    
    df %>% 
        mutate(lib_id = str_extract(lib_id, "lib[0-9]+"))
}

# clean/relabel donor ID
clean_donor_ids <- function(df) {
    df %>% 
        mutate(donor_id = tolower(donor_id))
}

# relabel timepoints
relabel_timepoints <- function(df) {
    df %>% 
        mutate(timepoint = str_replace(timepoint, " ", ""),
               timepoint = str_replace(timepoint, "InfusionProduct", "IP"),
               timepoint = str_replace(timepoint, "Day0", "t0"),
               timepoint = str_replace(timepoint, "Day(7|8|9|12)", "t1"),
               timepoint = str_replace(timepoint, "Day(26|28|29|33)", "t2"))
}

# convert transcript names to segment names
relabel_transcripts <- function(df, xcript_dat) {
    df_names <- names(df)
    for (i in 1:length(df_names)) {
        if (df_names[i] %in% xcript_dat$seqnames) {
            xcript_row <- which(xcript_dat$seqnames %in% df_names[i])
            new_name <- xcript_dat$transcript_id[xcript_row]
            df_names[i] <- new_name
        }
    }
    names(df) <- df_names
    return(df)
}