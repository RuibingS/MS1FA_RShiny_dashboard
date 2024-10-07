#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List adducts_anno_fun(Rcpp::List PI_res,
                            Rcpp::DataFrame FT,
                            double rt_thr,
                            double mz_thr,
                            int ppm) {
  
  int nrow_FT = FT.nrow();
  int length_PI_res = PI_res.length();
  
  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mz"];
  Rcpp::NumericVector FT_rt = FT["rt"];
  
  List output(length_PI_res);
  
  for(int i = 0; i < length_PI_res ; i++){
    Rcpp::List res_i = PI_res[i];  
    // Extract individual components from the list
    Rcpp::CharacterVector PI_feature_name = res_i["Feature_name"];
    Rcpp::NumericVector PI_Feature_rt = res_i["Feature_rt"];
    Rcpp::CharacterVector PI_name = res_i["PI_name"];
    Rcpp::NumericVector PI_M = res_i["M"];
    Rcpp::NumericVector PI_adducts_mass = res_i["adducts_mass"];
    Rcpp::CharacterVector PI_adducts_names = res_i["adduct_name"];

    

    // output
    Rcpp::List df(8);
    df.names() = Rcpp::CharacterVector::create("PI_feature_name_vec", "PI_name", "M", "adduct_name", "feature_name", "feature_mz", "feature_rt", "ion_mass");
    
    CharacterVector PI_featureNames;
    CharacterVector PINames;
    NumericVector MValues;
    CharacterVector adducts_Names;
    CharacterVector adducts_featureNames;
    Rcpp::NumericVector adducts_fea_mz;
    
    Rcpp::NumericVector adducts_mass_vec;
    Rcpp::NumericVector adducts_fea_RT;
    
    for(int k = 0; k < PI_adducts_mass.length(); k++){
      for(int j = 0; j < nrow_FT; j++){
        
        float diff_mz = std::fabs(FT_mz[j] - PI_adducts_mass[k]);
        float diff_rt = std::fabs(FT_rt[j] - PI_Feature_rt[0]);
        
        if((diff_mz <= PI_adducts_mass[k] * ppm * 1e-6 && diff_rt <= rt_thr) || (diff_mz <= mz_thr && diff_rt <= rt_thr)) {
          PI_featureNames.push_back(PI_feature_name[0]);
          PINames.push_back(PI_name[0]);
          MValues.push_back(PI_M[0]);
          adducts_Names.push_back(PI_adducts_names[k]);
          adducts_mass_vec.push_back(PI_adducts_mass[k]);
          adducts_featureNames.push_back(feature_name[j]);
          adducts_fea_mz.push_back(FT_mz[j]);
          adducts_fea_RT.push_back(FT_rt[j]);
        }
      }
    }
    
    df["PI_feature_name_vec"] = PI_featureNames;
    df["PI_name"] = PINames;
    df["M"] = MValues;
    df["adduct_name"] = adducts_Names;
    df["ion_mass"] = adducts_mass_vec;
    df["feature_name"] = adducts_featureNames;
    df["feature_mz"] = adducts_fea_mz;
    df["feature_rt"] = adducts_fea_RT;
    output[i] = df;
  }
  
  return output;
}
