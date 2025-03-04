###########################################
# MZmine output -> FDR                    #
###########################################

source(here::here("R","helper functions.R"))
source(here::here("R","parse_library.R"))

# step 1. read in MZmine feature table and Si16 library
# FDR for MZmine

# MZmine_FT <- read.csv(here::here("Data","feature_table","StM16_MZmine_export_FT_MS1.csv"))

# read Si16 libraray 
# Si16_test <- parse_library_file_parallel_NIST(file_path = here::here("Data","metabolite_data","Si16.library"), ionPolarity = "pos",spectrum_type ="MS1")



# step 2.split the feature table to list by the 16 standards mz and rt values

MZmine_split_FT_fun <- function(feature_table,tolerance) {

  feature_table_with_rt <- feature_table %>%
    mutate(rt = row.retention.time * 60,
           metabolite_name = trimws(sub("\\: \\[.*", "", row.identity..all.IDs.))
           ) %>%
    relocate(rt, .after = row.m.z) %>%
    data.frame()



  duplicated_metabolites <- feature_table_with_rt$metabolite_name[duplicated(feature_table_with_rt$metabolite_name)]
  duplicated_metabolites_filter <- duplicated_metabolites[which(duplicated_metabolites!="")]


  feature_table_with_rt_sub <- feature_table_with_rt %>%
    dplyr::group_by(metabolite_name) %>%
    dplyr::mutate(
      metabolite_name = if_else(
        metabolite_name %in% duplicated_metabolites_filter,
        paste0(trimws(sub("\\: \\[.*", "", metabolite_name)), "_", row_number()),
        trimws(sub("\\: \\[.*", "", metabolite_name))
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(metabolite_name != "" & !is.na(metabolite_name)) %>%
    data.frame()


  filtered_list <- list()


  unique_metabolites <- unique(feature_table_with_rt_sub$metabolite_name)
  for (metabolite in unique_metabolites) {

    metabolite_subset <- feature_table_with_rt_sub[feature_table_with_rt_sub$metabolite_name == metabolite, ]

    if (nrow(metabolite_subset) == 0) {
      next
    }

    reference_rt <- min(metabolite_subset$rt, na.rm = TRUE)

    matched_rows <- feature_table_with_rt[abs(feature_table_with_rt$rt - reference_rt) <= tolerance &
                                            !is.na(feature_table_with_rt$rt), ]

    filtered_list[[metabolite]] <- matched_rows
  }

  return(filtered_list)
}




# test
# MZmine_FT_list <- MZmine_split_FT_fun(feature_table = MZmine_FT,  tolerance = 3)



# Step 3. each compound and the nearby features matching to Si16 library

fuzzy_match_MZmine <- function(comp_list, FT_list, tolerance = 0.01) {

  FT_comp_names <- sub("_\\d+$", "", names(FT_list))

  output_temp <- setNames(lapply(comp_list, function(comp) {

    comp_name <- comp$Name
    comp_rt <- as.numeric(comp$RetTime)

    feature_list_temp <- FT_list[which(tolower(FT_comp_names) %in% tolower(comp_name))]

    if (length(feature_list_temp) == 0) {
      return(NULL)  
    }

    feature_list_temp_rt <- sapply(feature_list_temp, function(x) mean(x$rt, na.rm = TRUE))

    rt_diff <- abs(feature_list_temp_rt - comp_rt)
    closest_index <- which.min(rt_diff)

    df <- feature_list_temp[[closest_index]]

    df$library_name_found <- FALSE
    df$library_peak_matched <- FALSE
    df$check_group_id <- FALSE

    current_name_found <- tolower(comp_name) %in% tolower(FT_comp_names)

    if (current_name_found) {

      df <- df[, !grepl("Peak.area", colnames(df)), drop = FALSE]

      df_center_ID <- unique(df[which(grepl(tolower(comp_name),
                                            tolower(df$row.identity..all.IDs.),
                                            fixed = TRUE)), "correlation.group.ID"])

      if (length(df_center_ID[is.na(df_center_ID)]) > 1) {
        warning("More than 1 group index found!")
      }


      df$check_group_id <- ifelse(
        is.na(df$correlation.group.ID), 
        FALSE,
        ifelse(
          df$metabolite_name == "",  
          df$correlation.group.ID %in% df_center_ID,
          ifelse(
            comp_name != "" & tolower(comp_name) != tolower(df$metabolite_name),  
            FALSE,
            df$correlation.group.ID %in% df_center_ID  
          )
        )
      )
      
      
      
      df$library_name_found <- TRUE

      df$library_peak_matched <- sapply(seq_len(nrow(df)), function(row_idx) {
        mz_value <- df$row.m.z[row_idx]
        any(abs(comp$Peaks$mz - mz_value) <= tolerance)
      })
    }

    return(df)
  }), sapply(comp_list, function(comp) comp$Name))

  output_temp <- Filter(Negate(is.null), output_temp)

  return(output_temp)
}



# test output
# MZmine_result_merge <- fuzzy_match_MZmine(comp_list = Si16_test ,FT_list = MZmine_FT_list, tolerance = 0.01)
 
 
######################################################################
# FDR                                                                #
######################################################################

 
library(dplyr)
 
# Helper function
calculate_metrics <- function(df, consistency_col) {
   
   if (!all(c("library_name_found", "library_peak_matched", consistency_col) %in% colnames(df))) {
     stop(paste("Missing required columns in data:", consistency_col))
   }
   
   true_positives <- df %>%
     filter(library_name_found & library_peak_matched & !!sym(consistency_col)) %>%
     nrow()
   
   false_positives <- df %>%
     filter(library_name_found & !library_peak_matched & !!sym(consistency_col)) %>%
     nrow()
   
   false_negatives <- df %>%
     filter(library_name_found & library_peak_matched & ! (!!sym(consistency_col))) %>%
     nrow()
   
   true_negatives <- df %>%
     filter(library_name_found & !library_peak_matched & ! (!!sym(consistency_col))) %>%
     nrow()
   
   
   FDR <- ifelse((true_positives + false_positives) > 0,
                 false_positives / (true_positives + false_positives),
                 NA)
   
   TPR <- ifelse((true_positives + false_negatives) > 0,
                 true_positives / (true_positives + false_negatives),
                 NA)
   
   return(list(
     Confusion_Matrix = data.frame(
       True_Positive = true_positives,
       False_Positive = false_positives,
       False_Negative = false_negatives,
       True_Negative = true_negatives
     ),
     FDR = FDR,
     TPR = TPR
   ))
 }
 
# Main function for confusion matrix 

conf_matrix_MZmine <- function(MZmine_result_list ) 
 {
   
   
   df <- plyr::ldply(MZmine_result_list, data.frame, .id = "source")

   
   required_cols <- c("library_name_found", "library_peak_matched", "check_group_id")
   missing_cols <- setdiff(required_cols, colnames(df))
   if (length(missing_cols) > 0) {
     stop(paste("Missing columns in input data:", paste(missing_cols, collapse = ", ")))
   }
   
   
   results <- list(
     Group = calculate_metrics(df, "check_group_id")
     
     
   )
   
   return(results)
 }
 
 
 
 
# conf_matrix_MZmine_results <- conf_matrix_MZmine(MZmine_result_list= MZmine_result_merge)
# conf_matrix_MZmine_results

# New
 # $Group
 # $Group$Confusion_Matrix
 # True_Positive False_Positive False_Negative True_Negative
 # 1           268             48             82            52
 # 
 # $Group$FDR
 # [1] 0.1518987
 # 
 # $Group$TPR
 # [1] 0.7657143
 
 
# old
  # $Group
 # $Group$Confusion_Matrix
 # True_Positive False_Positive False_Negative True_Negative
 # 1           257             48             73            78
 # 
 # $Group$FDR
 # [1] 0.157377
 # 
 # $Group$TPR
 # [1] 0.7787879
 
#  $Group
#  $Group$Confusion_Matrix
#  True_Positive False_Positive False_Negative True_Negative
#  1           232             48             70           103
#  
#  $Group$FDR
#  [1] 0.1714286
#  
#  $Group$TPR
#  [1] 0.7682119
 























































































































