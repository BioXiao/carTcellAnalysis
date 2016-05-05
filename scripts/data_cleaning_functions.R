
# clean column names of data frame
clean_headers <- function(df) {
    df_names <- names(df)
    df_names <- str_to_lower(df_names)
    df_names <- str_replace_all(df_names, "( )+", "_")
    df_names <- str_replace_all(df_names, "[^[:alnum:]_]", "_")
    df_names <- str_replace_all(df_names, "_+", "_")
    df_names <- str_replace_all(df_names, "_$", "")
    
    names(df) <- df_names
    return(df)
}

# recursively deduplicate names in a vector
dedup_names <- function(x, n = 1) {
    n <- n + 1
    is_dup_name <- duplicated(x)
    x[is_dup_name] <- x[is_dup_name] %>% 
        map_chr(function(s) {
            ifelse(str_detect(s, str_c("_", n - 1)),
                   str_replace(s, str_c("_", n - 1), str_c("_", n)),
                   str_c(s, n, sep = "_"))
        })
    if(any(duplicated(x))) {
        dedup_names(x, n)
    } else {
        return(x)
    }
}

# clean up duplicated headers
clean_dup_names <- function(df) {
    df_names <- names(df)
    names(df) <- dedup_names(df_names)
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
        mutate(timepoint_short = str_replace(timepoint, " ", ""),
               timepoint_short = str_replace(timepoint_short, 
                                             "InfusionProduct", "IP"),
               timepoint_short = str_replace(timepoint_short, 
                                             "Day0", "t0"),
               timepoint_short = str_replace(timepoint_short, 
                                             "Day(7|8|9|12)", "t1"),
               timepoint_short = str_replace(timepoint_short, 
                                             "Day(26|28|29|33)", "t2"))
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