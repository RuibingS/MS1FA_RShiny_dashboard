#include <Rcpp.h>
#include <cmath> 
using namespace Rcpp;

// function to round a  number 
double roundToDecimalPlaces(double value, int places) {
  double multiplier = std::pow(10.0, places);
  return std::round(value * multiplier) / multiplier;
}

// [[Rcpp::export]]
DataFrame left_join_and_mutate_fun(DataFrame FT_df, 
                                   DataFrame long_df, double rt_thr) {
  CharacterVector fea_name = FT_df["feature_name"];
  NumericVector mz = FT_df["mz"];
  NumericVector rt = FT_df["rt"];
  
  CharacterVector FTname1 = long_df["Var1"];
  CharacterVector FTname2 = long_df["Var2"];
  NumericVector Freq_long = long_df["Freq"];
  
  int n = FTname1.size();

  NumericVector mz_x(n, NA_REAL), mz_y(n, NA_REAL);
  NumericVector rt_x(n, NA_REAL), rt_y(n, NA_REAL);
  NumericVector Freq_x(n, NA_REAL);
  NumericVector rt_diff(n), mz_diff(n);
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < fea_name.size(); ++j) {
      if (FTname1[i] == fea_name[j]) {
        mz_x[i] = mz[j];
        rt_x[i] = rt[j]; 
        Freq_x[i] = Freq_long[i];
      } 
      if (FTname2[i] == fea_name[j]) {
        mz_y[i] = mz[j];
        rt_y[i] = rt[j];
      }
    } 
    
    
    if (!NumericVector::is_na(mz_x[i])) mz_x[i] = roundToDecimalPlaces(mz_x[i], 5);
    if (!NumericVector::is_na(mz_y[i])) mz_y[i] = roundToDecimalPlaces(mz_y[i], 5);
    if (!NumericVector::is_na(Freq_x[i])) Freq_x[i] = roundToDecimalPlaces(Freq_x[i], 3);
    

    rt_diff[i] = std::abs(rt_x[i] - rt_y[i]);
    mz_diff[i] = std::abs(mz_x[i] - mz_y[i]);
  } 

  IntegerVector idx = seq(0, n - 1);
  idx = idx[rt_diff <= rt_thr];
  
  return DataFrame::create(_["Var1"] = FTname1[idx], _["Var2"] = FTname2[idx], _["Freq"] = Freq_x[idx], 
                           _["mz_x"] = mz_x[idx], _["mz_y"] = mz_y[idx], _["rt_x"] = rt_x[idx], 
                             _["rt_y"] = rt_y[idx], _["rt_diff"] = rt_diff[idx], _["mz_diff"] = mz_diff[idx],
                               _["stringsAsFactors"] = false);
} 

// [[Rcpp::export]]
DataFrame left_join_and_mutate_NL_fun(DataFrame FT_df, 
                                       DataFrame long_df, 
                                       double rt_thr) {
  CharacterVector fea_name = FT_df["feature_name"];
  NumericVector mz = FT_df["mz"];
  NumericVector rt = FT_df["rt"];
  
  CharacterVector FTname1=long_df["Var1"];
  CharacterVector FTname2=long_df["Var2"];
  
  

  int n = FTname1.size();
  
  // initialize result vectors
  NumericVector mz_x(n, NA_REAL), mz_y(n, NA_REAL);
  NumericVector rt_x(n, NA_REAL), rt_y(n, NA_REAL);
  NumericVector rt_diff(n), mz_diff(n);
  
  // Matching and calculating differences
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < fea_name.size(); ++j) {
      if (FTname1[i] == fea_name[j]) {
        mz_x[i] = mz[j];
        rt_x[i] = rt[j]; // Assuming you have a way to align these properly
       
      }
      if (FTname2[i] == fea_name[j]) {
        mz_y[i] = mz[j];
        rt_y[i] = rt[j];
      }
    } 
    
    mz_x[i] = roundToDecimalPlaces(mz_x[i], 5);
    mz_y[i] = roundToDecimalPlaces(mz_y[i], 5);
   
    
    rt_diff[i] = std::abs(rt_x[i] - rt_y[i]);
    mz_diff[i] = std::abs(mz_x[i] - mz_y[i]);
  }
  
  // filtering based on rt_diff
  IntegerVector idx = seq(0, n-1);
  idx = idx[rt_diff <= rt_thr];
  
  return DataFrame::create(_["Var1"] = FTname1[idx], _["Var2"] = FTname2[idx],
                           _["mz_x"] = mz_x[idx],_["mz_y"] = mz_y[idx],
                          _["rt_x"] = rt_x[idx], _["rt_y"] = rt_y[idx],
                          _["rt_diff"] = rt_diff[idx],_["mz_diff"] = mz_diff[idx], _["stringsAsFactors"] = false);
} 


// [[Rcpp::export]]

DataFrame left_join_MS2PI_fun(DataFrame FT_df, 
                                      DataFrame MS2PI_df
                                      ) {
  
  CharacterVector fea_name = FT_df["feature_name"];

  NumericVector median_intensity =FT_df["median_value"];
  
  CharacterVector MS2_FTname=MS2PI_df["feature_name"];
  CharacterVector PI_FTname=MS2PI_df["Precursor_matched_FT_name"];
  
  // Assuming rt and Freq vectors align with FTname1 and FTname2 by their orders
  int n = MS2_FTname.size();
  
  // Initialize result vectors

  NumericVector median_value_x(n, NA_REAL), median_value_y(n, NA_REAL);
  NumericVector median_value_diff(n, NA_REAL);
  LogicalVector matchFound(n, false);
  
  for (int i = 0; i < n; ++i) {
    
    for (int j = 0; j < fea_name.size(); ++j) {
      if (MS2_FTname[i] == fea_name[j]) {
        median_value_x[i] = median_intensity[j];
      } 
      if (PI_FTname[i] == fea_name[j]) {
        median_value_y[i] = median_intensity[j];
      } 
    }  
    // Calculate median_value_diff where possible
    if (!NumericVector::is_na(median_value_x[i]) && !NumericVector::is_na(median_value_y[i])) {
       double diff = median_value_y[i] - median_value_x[i];
      // if (diff >= 0) {
        median_value_diff[i] = diff;
        matchFound[i] = true; // Indicate a successful match with a positive difference
      // }
    }
  } 

    
  return DataFrame::create(_["feature_name"] = MS2_FTname[matchFound],
                           _["Precursor_matched_FT_name"] = PI_FTname[matchFound],
                           _["median_value_y"] = median_value_y[matchFound],
                           _["median_value_x"] = median_value_x[matchFound],
                           _["stringsAsFactors"] = false);
  
}





// [[Rcpp::export]]
DataFrame left_join_MS2PI_MZmine3_fun(DataFrame FT_df, DataFrame MS2PI_df) {
  
  CharacterVector fea_name = FT_df["feature_name"];
  NumericVector median_intensity = FT_df["median_value"];
  
  CharacterVector MS2_FTname = MS2PI_df["ms2_feature_name"];
  CharacterVector PI_FTname = MS2PI_df["PI_feature_name"];
  
  int n = MS2_FTname.size();
  
  NumericVector median_value_x(n, NA_REAL), median_value_y(n, NA_REAL);
  NumericVector median_value_diff(n, NA_REAL);
  LogicalVector matchFound(n, false);
  
  for (int i = 0; i < n; ++i) {
   
    for (int j = 0; j < fea_name.size(); ++j) {
      if (MS2_FTname[i] == fea_name[j]) {
        median_value_x[i] = median_intensity[j];
      }
      if (PI_FTname[i] == fea_name[j]) {
        median_value_y[i] = median_intensity[j];
      } 
    } 
    // Calculate median_value_diff where possible
    if (!NumericVector::is_na(median_value_x[i]) && !NumericVector::is_na(median_value_y[i])) {
      double diff = median_value_y[i] - median_value_x[i];
      if (diff >=0) {
        median_value_diff[i] = diff;
        matchFound[i] = true; // Indicate a successful match with a positive difference
      }
    }
  } 
  
  
  // Return a data frame with the filtered values
  return DataFrame::create(Named("ms2_feature_name") = MS2_FTname[matchFound],
                           Named("PI_feature_name") = PI_FTname[matchFound],
                           Named("median_value_x") = median_value_x[matchFound],
                           Named("median_value_y") = median_value_y[matchFound],
                           Named("stringsAsFactors") = false);
}

