#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector ISFAnno_fun(DataFrame mzMatch_df, DataFrame FT) {
  int n = mzMatch_df.nrows();
  int n_FT = FT.nrows();
  
  CharacterVector ms2_FTname = mzMatch_df["feature_name"];
  CharacterVector Precursor_FTname = mzMatch_df["Precursor_matched_FT_name"];
  CharacterVector FT_feature_name = FT["feature_name"];
  CharacterVector FT_ISF_anno(n_FT, "");
  
  //  check precursors in FT
  for (int j = 0; j < n_FT; j++) {
    for (int i = 0; i < n; i++) {
      if (FT_feature_name[j] == Precursor_FTname[i] ) { //&& ms2_FTname[i] == Precursor_FTname[i] ) {
        FT_ISF_anno[j] = "Precursor";
         break; 
      }
    }
  }
  
  // check FT features are both precursors and have MS2 matches
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n_FT; j++) {
      if (FT_feature_name[j] == ms2_FTname[i]) {
        std::string annotation = as<std::string>(Precursor_FTname[i]) + "_MS2 match";
        if (!FT_ISF_anno[j].empty() ) { //&& FT_ISF_anno[j] != "Precursor"
          FT_ISF_anno[j] = Rcpp::String(as<std::string>(FT_ISF_anno[j]) + "; " + annotation);
        } else { 
          FT_ISF_anno[j] = annotation;
        } 
        break;
      }
    }
  } 
  
  return FT_ISF_anno;
} 


  

    