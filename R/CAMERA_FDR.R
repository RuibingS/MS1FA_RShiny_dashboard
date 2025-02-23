# ###########################################
# #  XCMS output -> CAMERA -> FDR           #
# ###########################################

###################################################

source(here::here("R","helper functions.R"))
source(here::here("R","parse_library.R"))

# read in-house libraray 
# Si16_test <- parse_library_file_parallel_NIST(file_path = here::here("Data","metabolite_data","Si16_library.library"), ionPolarity = "pos",spectrum_type ="MS1")

# 
# # # read CAMERA output
# 
# 
# CAMERA_ft <- featureTable.import.fun(here::here("Data","feature_table","CAMERA_annotated_FT.csv"))



########################################################################################################
########################################################################################################

# Si11 library to find all features 

Si11_lib_FT_list <- function(comp_list, FT, rt_tolerance = 3, mz_tolerance = 0.01) {
  
  Si11_target_list <- metabolie.data.import.fun(here::here("Data","metabolite_data","new16Mix_targetlist.csv"), IonPolarity = "pos")
  # 
  if(!any("RT"==colnames(Si11_target_list))){
    Si11_target_list$RT <- rep(NA_real_, nrow(Si11_target_list))
  }
  
  PI_res <- PImatch_fun(FT = CAMERA_ft,Comp_data = Si11_target_list,ppm = 5,
                        PIon = c("[M+H]+","[M+Na]+"),diff_mz_thr = 0.002, diff_rt_thr= 10)
  
  PI_res_sub <- PI_res[which(sapply(PI_res, function(x) length(x$Feature_name)>0))]
  
  PI_match.df <- plyr::ldply(PI_res_sub, data.frame)
 
  MH.df_sub <- PI_match.df %>%
    dplyr::mutate(Comp_PI_name = paste(Comp_name, PI_name, sep = " ")) %>%
    dplyr::select(Feature_name,Comp_PI_name)
  
  
  # # self defined function in helper.r: convert long to wide format
  MH.df.output<-long_to_wide.fun(df=MH.df_sub,col_name1 =Feature_name,col_name2=Comp_PI_name )
  
  
  FT_anno <- FT %>% left_join(., MH.df.output, by = "feature_name") %>%   mutate(
               metabolite_name = trimws(sub("\\[.*", "", feature_annotation))) %>%data.frame()

  FT_comp_names <- sub("\\s*\\[.*$", "", FT_anno$feature_annotation)
  

  output_temp <- setNames(lapply(comp_list, function(comp) {
    
    comp_name <- comp$Name
    comp_rt <- as.numeric(comp$RetTime)

    matching_indices <- which(tolower(FT_comp_names) %in% tolower(comp_name))

    if (length(matching_indices) == 0) {
      return(NULL)
    }
    
    feature_df_temp <- FT_anno[matching_indices, ]
    
    rt_diff <- abs(feature_df_temp$rt - comp_rt)
    closest_index <- which.min(rt_diff)
    

    best_match <- feature_df_temp[closest_index, , drop = FALSE]
    
    FT_df_subset <- FT_anno[which(abs(FT_anno$rt - best_match$rt) <= rt_tolerance | 
                                    FT_anno$pcgroup %in% best_match$pcgroup), ]
    # 
    FT_df_subset$library_name_found <- FALSE
    FT_df_subset$library_peak_matched <- FALSE
    FT_df_subset$check_group_id <- FALSE
    
    current_name_found <- tolower(comp_name) %in% tolower(FT_comp_names)
    
    if (current_name_found) {
      
      df_center_ID <- unique(FT_df_subset[which(grepl(tolower(comp_name),
                                            tolower(FT_df_subset$feature_annotation),
                                            fixed = TRUE)), "pcgroup"])
   
      
      if (length(df_center_ID[is.na(df_center_ID)]) > 1) {
        warning("More than 1 group index found!")
      }
      
      
      FT_df_subset$check_group_id <- ifelse(
        is.na(FT_df_subset$pcgroup),
        FALSE,
        ifelse(
          is.na(FT_df_subset$metabolite_name),
          FT_df_subset$pcgroup %in% df_center_ID,
          ifelse(
            comp_name != "" & tolower(comp_name) != tolower(FT_df_subset$metabolite_name),
            FALSE,
            FT_df_subset$pcgroup %in% df_center_ID
          )
        )
      )
      

      FT_df_subset <- FT_df_subset[, !grepl("Peak.area", colnames(FT_df_subset)), drop = FALSE]

      FT_df_subset$library_name_found <- TRUE
      
      FT_df_subset$library_peak_matched <- sapply(FT_df_subset$mz, function(mz_value) {
        any(abs(comp$Peaks$mz - mz_value) <= mz_tolerance)
      })
    }
    
    return(FT_df_subset)
  }), sapply(comp_list, function(comp) comp$Name))
  output_temp <- Filter(Negate(is.null), output_temp)
  
  return(output_temp)
}

Si11_lib_FT_list_test <- Si11_lib_FT_list(comp_list = Si16_test, FT = CAMERA_ft, rt_tolerance = 3)


 
# ######################################################################
# # FDR                                                                #
# ######################################################################
# 
# 
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

# Main function
conf_matrix_CAMERA <- function(CAMERA_result_list) {
  
  df <- plyr::ldply(CAMERA_result_list, data.frame, .id = "source")
  
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



# CAMERA_results <- conf_matrix_CAMERA(CAMERA_result_list= Si11_lib_FT_list_test )
# CAMERA_results


 
# $Group
# $Group$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           431             79             30            14
# 
# $Group$FDR
# [1] 0.154902
# 
# $Group$TPR
# [1] 0.9349241








