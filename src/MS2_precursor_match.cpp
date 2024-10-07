#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List MS2match(Rcpp::S4 MS_OBJ,
                    Rcpp::DataFrame FT, 
                    double mz_diff_precursor,
                    double rt_thr_precursor,
                    double ppm)
{
  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mz"];
  Rcpp::NumericVector FT_RT = FT["rt"];
  
  int nrow_FT = FT.nrow();
  // the "ProtGenerics" package namespace
  Environment MSnbaseEnv = Environment::namespace_env("MSnbase");
   
  // the "spectra" method from the "ProtGenerics" package
  Function extractSpectraFun = MSnbaseEnv["extractSpectraData"];
  S4 spctr = extractSpectraFun(MS_OBJ);
   
  List reslist = spctr.slot("listData");
  NumericVector precursorMz = reslist["precursorMZ"];
  NumericVector precursorInt = reslist["precursorIntensity"];
  NumericVector rtime = reslist["rtime"];
  List mz = reslist["mz"];
  int nrow_Precursor = precursorMz.length();
  NumericVector precursor_mz_vec(nrow_Precursor);  // Initialize with the number of rows expected
  
  
  List output(nrow_FT);
  for (int i = 0; i < nrow_FT; i++) {
    Rcpp::List df(7);
    df.names() = Rcpp::CharacterVector::create("Feature_name", "Feature_mz", "Feature_RT", 
             "Precursor_mz","Precursor_intensity", "mz_Vals", "Precursor_rt");
    
    CharacterVector feaNames_vec;
    NumericVector feaRT_vec;
    NumericVector feature_mz_vec;
    NumericVector precursor_mz_vec;
    NumericVector precursor_intensity_vec;
    List mz_vals_list; // Changed to a list to store the mz_values vectors
    NumericVector Precursor_rt_vec;
    
    for (int j = 0; j < nrow_Precursor; j++) {
      double precursor_mz = precursorMz[j];
      double Precursor_rt = rtime[j];
      double precursor_intensity = precursorInt[j];
      List mz_values = mz[j];
      
      float diff_mz = std::fabs(FT_mz[i] - precursor_mz);
      float diff_rt = std::fabs(FT_RT[i] - Precursor_rt);
      
      if ((diff_mz <= precursor_mz * ppm * 1e-6 || diff_mz <= mz_diff_precursor)&& diff_rt <= rt_thr_precursor ) {
        feaNames_vec.push_back(feature_name[i]);
        feaRT_vec.push_back(FT_RT[i]);
        feature_mz_vec.push_back(FT_mz[i]);
        precursor_mz_vec.push_back(precursor_mz);
        precursor_intensity_vec.push_back(precursor_intensity);
        Precursor_rt_vec.push_back(Precursor_rt);
        mz_vals_list.push_back(mz_values); // Store the mz_values vector in the list
      }
    }
    
    df["Feature_name"] = feaNames_vec;
    df["Feature_mz"] = feature_mz_vec;
    df["Feature_RT"] = feaRT_vec;
    df["Precursor_mz"] = precursor_mz_vec;
    df["Precursor_rt"] = Precursor_rt_vec;
    df["Precursor_intensity"] = precursor_intensity_vec;
    df["mz_Vals"] = mz_vals_list; // Assign the list to mz_Vals
    output[i] = df;
  }
  
  return output;
} 
