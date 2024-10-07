#include <Rcpp.h>
using namespace Rcpp;


// Function to check if an object is a list
bool isList(SEXP x) {
  return Rf_isVectorList(x);
}

// [[Rcpp::export]]
List ISFMZmine_fun(Rcpp::List MZmine_list,
                   DataFrame FT,
                   double rt_thr_MS2,
                   double ppm,
                   double ms2_mz_diff
) {
  
  int nrow_FT = FT.nrows();
  CharacterVector FT_FTname = FT["feature_name"];
  NumericVector FT_rtime = FT["rt"];
  NumericVector FT_mz = FT["mz"];
  
  List output(MZmine_list.length());
  
  for (int i = 0; i < MZmine_list.length(); i++) {
    
    if (isList(MZmine_list[i])) {
      
      Rcpp::List inner_list = MZmine_list[i];
      CharacterVector PI_FTname = inner_list["feature_name"];
      NumericVector PI_rtime = inner_list["rt"];
      NumericVector PI_mz = inner_list["PI_mass"];
      Rcpp::List ms2df = inner_list["MS2Df"];
      NumericVector ms2mz = ms2df["mz"];
      
      Rcpp::List df(7);
      df.names() = Rcpp::CharacterVector::create("PI_feature_name","PI_feature_mz","PI_feature_rt","ms2_feature_name", "ms2_feature_mz","ms2_mz", "ms2_feature_rt");
      
      CharacterVector PI_feature_name_temp;
      NumericVector PI_feature_mz_temp;
      NumericVector PI_feature_rt_temp;
      CharacterVector ms2_feature_name_temp;
      NumericVector FTfeature_mz_temp;
      NumericVector ms2_feature_mz_temp;
      NumericVector ms2_feature_rt_temp;

      
      for (int j = 0; j < nrow_FT; j++) {
        
        for (int k = 0; k < ms2mz.length(); k++) {
          
          double mz_diff = std::abs(ms2mz[k] - FT_mz[j]);
          double rt_diff = std::abs(PI_rtime[0] - FT_rtime[j]); // Check if this index is correct
          
          if ((mz_diff <= ms2_mz_diff && rt_diff <= rt_thr_MS2) || (mz_diff <= ms2mz[k] * ppm * 1e-6 && rt_diff <= rt_thr_MS2)) {
            PI_feature_name_temp.push_back(PI_FTname[0]);
            PI_feature_mz_temp.push_back(PI_mz[0]);
            PI_feature_rt_temp.push_back(PI_rtime[0]);
            FTfeature_mz_temp.push_back(FT_mz[j]);
            ms2_feature_name_temp.push_back(FT_FTname[j]);
            ms2_feature_mz_temp.push_back(ms2mz[k]);
            ms2_feature_rt_temp.push_back(FT_rtime[j]);
          }
        }
      }
      
      df["PI_feature_name"] = PI_feature_name_temp;
      df["PI_feature_mz"] = PI_feature_mz_temp;
      df["PI_feature_rt"] = PI_feature_rt_temp;
      df["ms2_feature_mz"] = FTfeature_mz_temp;
      df["ms2_feature_name"] = ms2_feature_name_temp;
      df["ms2_mz"] = ms2_feature_mz_temp;
      df["ms2_feature_rt"] = ms2_feature_rt_temp;
      
      output[i] = df;
    }
  }
  return output;
}

// [[Rcpp::export]]
CharacterVector ISFMZmine_assign_fun(DataFrame PI_MS2_df,
                                     DataFrame FT) {
  
  int n = PI_MS2_df.nrows();
  int n_FT = FT.nrows();
  
  // MS2 match table
  CharacterVector ms2_FTname = PI_MS2_df["ms2_feature_name"];
  CharacterVector Precursor_FTname = PI_MS2_df["PI_feature_name"];
  
  // feature table
  CharacterVector FT_feature_name = FT["feature_name"];
  CharacterVector FT_ISF_anno(n_FT, "");
  
  //  check precursors in FT
  for (int j = 0; j < n_FT; j++) {
    for (int i = 0; i < n; i++) {
      if (FT_feature_name[j] == Precursor_FTname[i]) {
        FT_ISF_anno[j] = "Precursor";
        break; 
      }
    }
  } 
  
  // assign MS2 matches both precursors and have MS2 matches
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n_FT; j++) {
      if (FT_feature_name[j] == ms2_FTname[i]) {
        if (!FT_ISF_anno[j].empty()) {
          FT_ISF_anno[j] = Rcpp::String(as<std::string>(FT_ISF_anno[j]) + "; " + as<std::string>(Precursor_FTname[i]) + "_MS2 match");
        } else { 
          FT_ISF_anno[j] = as<std::string>(Precursor_FTname[i]) + "_MS2 match";
        } 
        break; 
      }
    }
  }
   
  return FT_ISF_anno;
} 