
# format for feature table output with feature annotation
# helper function: stretch long to wide format then collapse columns in to one cell
long_to_wide.fun<-function(df, col_name1, col_name2){
  
  
  output.df<-df%>%
    dplyr::select({{col_name1}},{{col_name2}}) %>%
    dplyr::group_by({{col_name1}}) %>%
    distinct() %>%
    dplyr::mutate(Var= paste0("annotation_", dplyr::row_number())) %>%
    tidyr::spread(Var,{{col_name2}} ,fill=NA_character_) %>%
    data.frame()
  
  output<-tidyr::unite(dplyr::distinct(output.df), col='fea_annotation', contains("annotation"), sep='; ',na.rm = TRUE)
  
  colnames(output)<-c("feature_name","feature_annotation")
  
  return(output)
}

#########################################################################
##############################################################################################
# helper function: intensity values row wise scale
df_scale.fun<-function(featureTable, start_col,end_col){
  
  start_col_ind<-which(colnames(featureTable)==start_col)
  end_col_ind<-which(colnames(featureTable)==end_col)
  
  df_row.scale<-t(apply(featureTable[, start_col_ind:end_col_ind], 1, scale))
  colnames(df_row.scale)<-colnames(featureTable[, start_col_ind:end_col_ind])
  rownames(df_row.scale)<-featureTable$feature_name
  t_df_row.scale<-t(df_row.scale)
  colnames(t_df_row.scale)
  return(t_df_row.scale)
}#


##############################################################################
## helper function
## input: data frame which has column names: feature_name and feature_annotation
## if one feature_name has multiple matches, then paste into one column.



feature_anno_rename_fun<-function(input_df){
  
  output_wide<-input_df %>%
    select(feature_name,feature_annotation) %>%
    group_by(feature_name) %>%
    distinct() %>%
    dplyr::mutate(Var= paste0("annotation_", dplyr::row_number())) %>%
    tidyr::spread(Var,feature_annotation ,fill=NA_character_) %>%
    data.frame()
  
  output<-tidyr::unite(output_wide, col='feature_annotation', contains("annotation"), sep='; ',na.rm = TRUE,)
  
  
  return(output)
  
}
########################################################################################
# calculate the adducts ion mass


adducts_mass_function<-function(M,adduct_formula){
  
  # remove space
  if (any(grepl(" ",adduct_formula))) {
    ion_mass<-gsub(" ", "",adduct_formula)
  }
  else{
    ion_mass<-adduct_formula
  }
  
  # Define the vector of expressions
  expr_vector <-ion_mass
  
  # Function to evaluate expressions for a given value of m
  evaluate_expr <- function(expr, M) {
    expr_with_multiplication <- gsub("(\\d+)M", "\\1 * M", expr)
    
    expr <- parse_expr(expr_with_multiplication)
    eval_tidy(expr, list(M = M))
  }
  
  result_vector <- sapply(expr_vector, function(expr) evaluate_expr(expr, M))
  
  
  return(as.vector(result_vector))
  
}

# test adducts formula


########################################################################################
# check if a compound feature is another compound's ISF feature 
#######################################################################################
library(dplyr)
library(stringr)
library(purrr)

check_ISF_fun <- function(df) {
  df %>%
    mutate(
      isf_refs = map(ISF_anno, ~ str_extract_all(.x, "FT\\d+_MS2 match")[[1]] %>%
                       str_remove("_MS2 match"))
    ) %>%
    mutate(
      metabolite_refs = map(isf_refs, ~ df %>%
                              filter(feature_name %in% .x) %>%
                              pull(metabolite_annotation) %>%
                              na.omit() %>%
                              str_split(";\\s*") %>%
                              unlist() %>%
                              str_extract("^[^\\[]*") %>%
                              str_trim()
      )
    ) %>%
    mutate(
      annotation_compounds = map(metabolite_annotation, ~ str_split(.x, ";\\s*") %>%
                                   unlist() %>%
                                   str_extract("^[^\\[]*") %>%
                                   str_trim()
      )
    ) %>%
    mutate(
      metabolite_annotation = if_else(
        map2_lgl(metabolite_refs, annotation_compounds, ~ any(!.x %in% .y)),
        NA_character_,
        metabolite_annotation
      ),
      adducts_anno = if_else(
        map2_lgl(metabolite_refs, annotation_compounds, ~ any(!.x %in% .y)),
        NA_character_,
        adducts_anno
      )
    ) %>%
    select(-isf_refs, -metabolite_refs, -annotation_compounds)  # Remove temporary columns
}


################################
# helper function

check_consistency_fun <- function(data) {
  
  data$metabolite_name <- trimws(sub("\\s*\\[M.*", "", data$metabolite_annotation))
  
  data <- data %>%
    mutate(
      check_group_consistency = FALSE,
      check_cor_group_consistency = FALSE,
      check_both_group_consistency = FALSE
    )
  
  unique_metabolites <- unique(data$metabolite_name[!is.na(data$metabolite_name)])
  
  for (metabolite in unique_metabolites) {
    
    temp_data <- data %>% filter(metabolite_name == metabolite)
    
    unique_groups <- unique(temp_data$group[!is.na(temp_data$group)])
    
    unique_cor_groups <- unique(temp_data$cor_group[!is.na(temp_data$cor_group)])
    
    if (length(unique_groups) > 1) {
      stop(paste("Error: More than one unique group associated with metabolite:", metabolite))
    }
    
    # Update consistency checks
    reference_group <- unique_groups
    reference_cor_group <- unique_cor_groups
    
    data <- data %>%
      mutate(
        check_group_consistency = if_else(
          group %in% reference_group,
          TRUE,
          check_group_consistency
        ),
        check_cor_group_consistency = if_else(
          cor_group %in% reference_cor_group,
          TRUE,
          check_cor_group_consistency
        ),
        check_both_group_consistency = if_else(
          check_group_consistency & check_cor_group_consistency,
          TRUE,
          check_both_group_consistency
        ),
        check_Or_group_consistency = if_else(
          check_group_consistency | check_cor_group_consistency,
          TRUE,
          check_both_group_consistency
        )
      )
  }
  
  return(data)
}

