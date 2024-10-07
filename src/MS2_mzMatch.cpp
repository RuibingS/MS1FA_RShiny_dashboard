#include <Rcpp.h>
using namespace Rcpp;

// check if an object is a list
bool isList(SEXP x) {
  return Rf_isVectorList(x);
}

// [[Rcpp::export]]
Rcpp::List MS2_mzMatch(Rcpp::List input_list, 
                       Rcpp::DataFrame FT,
                       double mz_diff_MS2,
                       double rt_thr_MS2,
                       double ppm
) {
  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mz"];
  Rcpp::NumericVector FT_RT = FT["rt"];
  int nrow_FT = FT.nrow();
  int length_list = input_list.length();
  
  Rcpp::List results; // Create an empty list to store the results
  
  for (int i = 0; i < nrow_FT; i++) {
    for (int j = 0; j < length_list; j++) {
      // Check if the element is a list
      if (isList(input_list[j])) {
        
        Rcpp::List inner_list = input_list[j];
        Rcpp::List mz_Vals_list = inner_list["mz_Vals"];
        int length_mz_val = mz_Vals_list.length();
        
        // Check if the inner list is empty
        if (length_mz_val > 0) {
          Rcpp::CharacterVector Feature_name_vec = inner_list["Feature_name"];
          Rcpp::NumericVector Feature_RT_vec= inner_list["Feature_RT"];
          Rcpp::NumericVector Precursor_rt_vec = inner_list["Precursor_rt"];
          Rcpp::NumericVector Precursor_mz_vec = inner_list["Precursor_mz"];
          
          for (int k = 0; k < length_mz_val; k++) {
            Rcpp::List mz_vals_inner_list = Rcpp::as<Rcpp::List>(mz_Vals_list[k]);
            int mz_vals_inner_list_length = mz_vals_inner_list.length();
            
            for (int m = 0; m < mz_vals_inner_list_length; m++) {
              Rcpp::NumericVector mz_values = Rcpp::as<Rcpp::NumericVector>(mz_vals_inner_list[m]);
              int mz_values_length = mz_values.length();
              
              for(int v = 0; v < mz_values_length; v++){
                if (mz_values_length > 0 && (std::abs(mz_values[v] - FT_mz[i]) <= mz_diff_MS2 || std::abs(mz_values[v] - FT_mz[i]) <=  FT_mz[i] * ppm * 1e-6)&& std::abs(Feature_RT_vec[0]-FT_RT[i])<=rt_thr_MS2 ) { //
                  Rcpp::List result; // Create a list to store the result for this match
                  
                  result["feature_name"] = Rcpp::as<std::string>(feature_name[i]);
                  result["feature_mz"] =FT_mz[i] ;
                  result["feature_rt"] =FT_RT[i] ;
                  result["Precursor_matched_FT_name"] = Rcpp::as<std::string>(Feature_name_vec[0]);
                  result["Precursor_mz"] = Precursor_mz_vec[0];
                  result["Precursor_matched_FT_rt"] =Feature_RT_vec[0];
                  
                  result["Precursor_rt"] = Precursor_rt_vec[0];
                  result["found_MS2_mz"] = mz_values[v];
                  results.push_back(result); // Add the result to the list of results
                }
              }
            }
          }
        }
      }
    }
  }
  
  return results; // Return the list of results
}
