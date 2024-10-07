
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







