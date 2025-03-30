############################################
#   MS1FA output -> FDR                    #
############################################

source(here::here("R","helper functions.R"))
source(here::here("R","parse_library.R"))

# read in Si16 library 
# Si16_test <- parse_library_file_parallel_NIST(file_path = here::here("Data","metabolite_data","Si16.library"), ionPolarity = "pos",spectrum_type ="MS1")
 
#####################################################
# 1. MS1FA output FT by XCMS processing             #
#####################################################

# import MS1FA output feature table 


# MS1FA_xcms_FT <- read.csv(here::here("output","MS1FA_output_FT","MS1FA_Feature_table_output_XCMS.csv"))

# MS1FA_mzmine_FT <- read.csv(here::here("output","MS1FA_output_FT","MS1FA_Feature_table_output_MZmine.csv")) 

MS1FA_MS1match_FT <- read.csv(here::here("output","MS1FA_output_FT","Feature_table_output_MS1match.csv")) 



# confusion matrix function


library(dplyr)

conf_matrix_fun <- function(input_list )
  { 
  
  merge_df <- plyr::ldply(input_list, data.frame)
  

  calculate_metrics <- function(df, consistency_col) {
    
    true_positives <- merge_df %>% 
      filter(library_name_found & library_peak_matched & !!sym(consistency_col)) %>% 
      nrow()
    
    false_positives <- merge_df %>% 
      filter(!library_peak_matched & !!sym(consistency_col)) %>% 
      nrow()
    
    false_negatives <- merge_df %>% 
      filter(library_name_found & library_peak_matched & ! (!!sym(consistency_col))) %>% 
      nrow()
    
    true_negatives <- merge_df %>% 
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
  
  results <- list(
    check_group_id = calculate_metrics(merge_df, "check_group_id"),
    OR = calculate_metrics(merge_df, "check_Or_group_consistency"),
    AND = calculate_metrics(merge_df, "check_both_group_consistency"),
    Group = calculate_metrics(merge_df, "check_group_consistency"),
    CorGroup = calculate_metrics(merge_df, "check_cor_group_consistency")
  )
  
  return(results)
}




MS1FA_FT_split <- function(comp_list, FT, rt_tolerance, mz_tolerance = 0.01) {
  
  # target list 
  Si16_target_list <- metabolie.data.import.fun(here::here("Data","metabolite_data","new16Mix_targetlist.csv"), IonPolarity = "pos")
  
  if(!any("RT"==colnames(Si16_target_list))){
    Si16_target_list$RT <- rep(NA_real_, nrow(Si16_target_list))
  }
  
  # PI match 
  PI_res <- PImatch_fun(FT,Comp_data = Si16_target_list,ppm = 5,
                        PIon = c("[M+H]+","[M+Na]+"),diff_mz_thr = 0.002, diff_rt_thr= 5)
  
  PI_res_sub <- PI_res[which(sapply(PI_res, function(x) length(x$Feature_name)>0))]
  
  PI_match.df <- plyr::ldply(PI_res_sub, data.frame)
  
  
  MH.df_sub <- PI_match.df %>%
    dplyr::mutate(Comp_PI_name = paste(Comp_name, PI_name, sep = " ")) %>%
    dplyr::select(Feature_name,Comp_PI_name)
  
  # helper function in helper.r: convert long to wide format
  MH.df.output<-long_to_wide.fun(df=MH.df_sub,col_name1 =Feature_name,col_name2=Comp_PI_name )
  
  FT_anno <- FT %>% 
    left_join(., MH.df.output, by = "feature_name") %>%   
    mutate(metabolite_name = trimws(sub("\\[.*", "", feature_annotation))) %>%
    data.frame()
  
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
    
    FT_df_subset <- FT_anno[which(abs(FT_anno$rt - best_match$rt) <= rt_tolerance), ]
    
    if (nrow(FT_df_subset) == 0) {
      return(NULL)
    }

    FT_df_subset$library_name_found <- FALSE
    FT_df_subset$library_peak_matched <- FALSE
    FT_df_subset$check_group_id <- FALSE
    
    current_name_found <- tolower(comp_name) %in% tolower(FT_comp_names)
    
    if (current_name_found) {
      
      df_center_ID <- unique(FT_df_subset$group[grepl(tolower(comp_name),
                                                      tolower(FT_df_subset$feature_annotation),
                                                      fixed = TRUE)])
      df_center_ID_cor <- unique(FT_df_subset$cor_group [grepl(tolower(comp_name),
                                                             tolower(FT_df_subset$feature_annotation),
                                                             fixed = TRUE)])
      

      if (length(df_center_ID) == 0) df_center_ID <- NA
      if (length(df_center_ID_cor) == 0) df_center_ID_cor <- NA
      
      if (sum(is.na(df_center_ID)) > 1) {
        warning("More than 1 group index found in df_center_ID!")
      }
      if (sum(is.na(df_center_ID_cor)) > 1) {
        warning("More than 1 group index found in df_center_ID_cor!")
      }
      
      if (nrow(FT_df_subset) > 0) {
        FT_df_subset$check_group_id <- ifelse(
          (!is.na(FT_df_subset$group) & FT_df_subset$group %in% df_center_ID) | 
            (!is.na(FT_df_subset$cor_group) & FT_df_subset$cor_group %in% df_center_ID_cor), 
          TRUE, 
          FALSE
        )
        FT_df_subset$center_name <- ifelse(
          (!is.na(FT_df_subset$group) & FT_df_subset$group %in% df_center_ID) | 
            (!is.na(FT_df_subset$cor_group) & FT_df_subset$cor_group %in% df_center_ID_cor), 
          tolower(comp_name), 
          FALSE
        )
        
      }
      
      FT_df_subset <- FT_df_subset[, !grepl("Peak.area", colnames(FT_df_subset)), drop = FALSE]
      
      FT_df_subset$library_name_found <- TRUE
      
      FT_df_subset$library_peak_matched <- sapply(FT_df_subset$mz, function(mz_value) {
        any(abs(comp$Peaks$mz - mz_value) <= mz_tolerance)
      })
    }
    
    FT_output <- check_consistency_fun(FT_df_subset)
    
    return(FT_output)
  }), sapply(comp_list, function(comp) comp$Name))
  
  output_temp <- Filter(Negate(is.null), output_temp)
  
  return(output_temp)
}

# test XCMS FT split
# Si16_lib_FT_list_MS1FA_rt3 <- MS1FA_FT_split(comp_list = Si16_test, FT = MS1FA_xcms_FT, rt_tolerance = 3)

Si16_lib_FT_list_MS1FA_rt3 <- MS1FA_FT_split(comp_list = Si16_test, FT = MS1FA_MS1match_FT, rt_tolerance = 3)

conf_matrix_MS1FA_XCMS <- conf_matrix_fun(Si16_lib_FT_list_MS1FA_rt3)
conf_matrix_MS1FA_XCMS


####################################################################################
# conf_matrix_MS1FA_XCMS <- conf_matrix_fun(Si16_lib_FT_list_MS1FA_rt3)
# conf_matrix_MS1FA_XCMS
# $check_group_id
# $check_group_id$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           437             44             24            53
# 
# $check_group_id$FDR
# [1] 0.09147609
# 
# $check_group_id$TPR
# [1] 0.9479393
# 
# 
# $OR
# $OR$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           437             44             24            53
# 
# $OR$FDR
# [1] 0.09147609
# 
# $OR$TPR
# [1] 0.9479393
# 
# 
# $AND
# $AND$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           313             10            148            87
# 
# $AND$FDR
# [1] 0.03095975
# 
# $AND$TPR
# [1] 0.6789588
# 
# 
# $Group
# $Group$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           375             34             86            63
# 
# $Group$FDR
# [1] 0.08312958
# 
# $Group$TPR
# [1] 0.813449
# 
# 
# $CorGroup
# $CorGroup$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           375             20             86            77
# 
# $CorGroup$FDR
# [1] 0.05063291
# 
# $CorGroup$TPR
# [1] 0.813449
# 







#####################################################
# 2. MS1FA output FT by MZmine processing           #
#####################################################

MS1FA_FT_MZmine_split <- function(comp_list, FT, rt_tolerance, mz_tolerance = 0.01) {
  
  Si16_target_list <- metabolie.data.import.fun(here::here("Data","metabolite_data","new16Mix_targetlist.csv"), IonPolarity = "pos")
  if(!any("RT"==colnames(Si16_target_list))){
    Si16_target_list$RT <- rep(NA_real_, nrow(Si16_target_list))
  }
  
  
  PI_res_MZmine <- PImatch_fun(FT,Comp_data = Si16_target_list,ppm = 5,
                               PIon = c("[M+H]+","[M+Na]+"),diff_mz_thr = 0.002, diff_rt_thr= 10)
  
  PI_res_MZmine_sub <- PI_res_MZmine[which(sapply(PI_res_MZmine, function(x) length(x$Feature_name)>0))]
  
  PI_match_MZmine.df <- plyr::ldply(PI_res_MZmine_sub, data.frame)
  
  MZmine_MH.df_sub <- PI_match_MZmine.df %>%
    dplyr::mutate(Comp_PI_name = paste(Comp_name, PI_name, sep = " ")) %>%
    dplyr::select(Feature_name,Comp_PI_name)
  
  
  # Helper function in helper.r: convert long to wide format
  MZmine_MH.df.output<-long_to_wide.fun(df=MZmine_MH.df_sub,col_name1 =Feature_name,col_name2=Comp_PI_name )
  
  
  FT_anno <- FT %>% 
    left_join(., MZmine_MH.df.output, by = "feature_name") %>%   
    mutate(metabolite_name = trimws(sub("\\[.*", "", feature_annotation))) %>%
    data.frame()
  
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
    
    FT_df_subset <- FT_anno[which(abs(FT_anno$rt - best_match$rt) <= rt_tolerance), ]
    
    if (nrow(FT_df_subset) == 0) {
      return(NULL)
    }
    
    FT_df_subset$library_name_found <- FALSE
    FT_df_subset$library_peak_matched <- FALSE
    FT_df_subset$check_group_id <- FALSE
    
    current_name_found <- tolower(comp_name) %in% tolower(FT_comp_names)
    
    if (current_name_found) {
      
      df_center_ID <- unique(FT_df_subset$group[grepl(tolower(comp_name),
                                                      tolower(FT_df_subset$feature_annotation),
                                                      fixed = TRUE)])
      df_center_ID_cor <- unique(FT_df_subset$cor_group [grepl(tolower(comp_name),
                                                               tolower(FT_df_subset$feature_annotation),
                                                               fixed = TRUE)])
      
      
      if (length(df_center_ID) == 0) df_center_ID <- NA
      if (length(df_center_ID_cor) == 0) df_center_ID_cor <- NA
  
      if (sum(is.na(df_center_ID)) > 1) {
        warning("More than 1 group index found in df_center_ID!")
      }
      if (sum(is.na(df_center_ID_cor)) > 1) {
        warning("More than 1 group index found in df_center_ID_cor!")
      }
      
      

      if (nrow(FT_df_subset) > 0) {
        FT_df_subset$check_group_id <- ifelse(
          (!is.na(FT_df_subset$group) & FT_df_subset$group %in% df_center_ID) | 
            (!is.na(FT_df_subset$cor_group) & FT_df_subset$cor_group %in% df_center_ID_cor), 
          TRUE, 
          FALSE
        )
      }
      
      FT_df_subset <- FT_df_subset[, !grepl("Peak.area", colnames(FT_df_subset)), drop = FALSE]
      
      FT_df_subset$library_name_found <- TRUE
      
      FT_df_subset$library_peak_matched <- sapply(FT_df_subset$mz, function(mz_value) {
        any(abs(comp$Peaks$mz - mz_value) <= mz_tolerance)
      })
    }
    
    FT_output <- check_consistency_fun(FT_df_subset)
    
    return(FT_output)
  }), sapply(comp_list, function(comp) comp$Name))
  
  output_temp <- Filter(Negate(is.null), output_temp)
  
  return(output_temp)
}


# test: split and match to in-house library

# Si16_lib_FT_list_MS1FA_mzmine <- MS1FA_FT_MZmine_split(comp_list = Si16_test, FT = MS1FA_mzmine_FT, rt_tolerance = 3)



# test ouotput: confusion matrix 
# conf_matrix_MS1FA_mzmine <- conf_matrix_fun(Si16_lib_FT_list_MS1FA_mzmine)
# 
# conf_matrix_MS1FA_mzmine

# 
# $check_group_id
# $check_group_id$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           316             42             34            58
# 
# $check_group_id$FDR
# [1] 0.1173184
# 
# $check_group_id$TPR
# [1] 0.9028571
# 
# 
# $OR
# $OR$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           317             43             33            57
# 
# $OR$FDR
# [1] 0.1194444
# 
# $OR$TPR
# [1] 0.9057143
# 
# 
# $AND
# $AND$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           210             14            140            86
# 
# $AND$FDR
# [1] 0.0625
# 
# $AND$TPR
# [1] 0.6
# 
# 
# $Group
# $Group$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           259             29             91            71
# 
# $Group$FDR
# [1] 0.1006944
# 
# $Group$TPR
# [1] 0.74
# 
# 
# $CorGroup
# $CorGroup$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           268             28             82            72
# 
# $CorGroup$FDR
# [1] 0.09459459
# 
# $CorGroup$TPR
# [1] 0.7657143


# Old
# output
# Si16_lib_FT_list_MS1FA_mzmine_res
# $check_group_id
# $check_group_id$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           308             42             28            57
# 
# $check_group_id$FDR
# [1] 0.12
# 
# $check_group_id$TPR
# [1] 0.9166667
# 
# 
# $OR
# $OR$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           308             42             28            57
# 
# $OR$FDR
# [1] 0.12
# 
# $OR$TPR
# [1] 0.9166667
# 
# 
# $AND
# $AND$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           203             13            133            86
# 
# $AND$FDR
# [1] 0.06018519
# 
# $AND$TPR
# [1] 0.6041667
# 
# 
# $Group
# $Group$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           250             28             86            71
# 
# $Group$FDR
# [1] 0.1007194
# 
# $Group$TPR
# [1] 0.7440476
# 
# 
# $CorGroup
# $CorGroup$Confusion_Matrix
# True_Positive False_Positive False_Negative True_Negative
# 1           261             27             75            72
# 
# $CorGroup$FDR
# [1] 0.09375
# 
# $CorGroup$TPR
# [1] 0.7767857
# 






