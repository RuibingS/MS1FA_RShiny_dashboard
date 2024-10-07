#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// "[M+H]+","[M+Na]+" PI match
Rcpp::List PImatch_fun(Rcpp::DataFrame FT,
                            Rcpp::DataFrame Comp_data,
                            Rcpp::CharacterVector PIon,
                            double diff_mz_thr,
                            double diff_rt_thr,
                            int ppm) {
  int nrow_FT = FT.nrow();
  int nrow_Comp = Comp_data.nrow();
  int length_PIon = PIon.length();

  Rcpp::CharacterVector feature_name = FT["feature_name"];
  Rcpp::NumericVector FT_mz = FT["mz"];
  Rcpp::NumericVector FT_rt = FT["rt"];


  Rcpp::CharacterVector Comp_name = Comp_data["Name"];
  Rcpp::NumericVector Comp_rt = Comp_data["RT"];
  Rcpp::NumericVector Comp_M = Comp_data["exactMass"];
  
  List output(nrow_FT);
  for(int i = 0; i < nrow_FT ; i++)
  {
    Rcpp::List df(6);
    df.names() = Rcpp::CharacterVector::create("Feature_name", "Feature_mz","Feature_rt","Comp_name","PI_name","M");

    CharacterVector feaNames;
    CharacterVector CompNames;
    CharacterVector PIname;
    Rcpp::NumericVector  fea_rt;
    Rcpp::NumericVector  fea_mz;
    Rcpp::NumericVector  M_vec;

    for(int p=0; p<length_PIon ; p++){
      String PIon_name = PIon[p];
      Rcpp::NumericVector comp_mz = Comp_data[PIon_name];


      for(int j = 0; j < nrow_Comp; j++){


        float diff_mz= std::fabs(FT_mz[i]-comp_mz[j]);
        
        if ((Comp_rt[j] != 0 || !NumericVector::is_na(Comp_rt[j]))&&diff_mz<= comp_mz[j]*ppm* 1e-6 && std::abs(FT_rt[i]-Comp_rt[j])<=diff_rt_thr||
             diff_mz<= diff_mz_thr && std::abs(FT_rt[i]-Comp_rt[j])<=diff_rt_thr) { 
          //{
          
          feaNames.push_back(feature_name[i]);
          fea_rt.push_back(FT_rt[i]);
          fea_mz.push_back(FT_mz[i]);
          CompNames.push_back(Comp_name[j]);
          PIname.push_back(PIon[p]);
          M_vec.push_back(Comp_M[j]);
          
          //}

        }
        
        else{
            if(diff_mz<= comp_mz[j]*ppm* 1e-6||diff_mz<= diff_mz_thr ) 
        {
          feaNames.push_back(feature_name[i]);
          fea_rt.push_back(FT_rt[i]);
          fea_mz.push_back(FT_mz[i]);
          CompNames.push_back(Comp_name[j]);
          PIname.push_back(PIon[p]);
          M_vec.push_back(Comp_M[j]);
          
        }
        }

      }

    }

    df["Feature_name"] = feaNames;
    df["Feature_mz"] = fea_mz;
    df["Feature_rt"] = fea_rt;
    df["Comp_name"] = CompNames;
    df["PI_name"] = PIname;
    df["M"] = M_vec;
    output[i] = df;
  }
  return output;

}

// 
// // [[Rcpp::export]]
// // NL match
// Rcpp::List findNL_fun(
//     Rcpp::DataFrame NL_data,
//     Rcpp::DataFrame adj_long,
//     double diff_mz_thr,
//     double diff_rt_thr,
//     int ppm
//     ) {
// 
// 
// 
//   int nrow_adj_long = adj_long.nrow();
//   int nrow_NL = NL_data.nrow();
// 
//   Rcpp::NumericVector mz_diff=adj_long["mz_diff"];
//   Rcpp::NumericVector rt_diff=adj_long["rt_diff"];
//   Rcpp::NumericVector NL_mz= NL_data["Accurate.Mass"];
//   Rcpp::CharacterVector NL_name = NL_data["Neutral.Loss"];
//   
// 
//   Rcpp::NumericVector mz1 = adj_long["mz_x"]; 
//   Rcpp::NumericVector mz2 = adj_long["mz_y"]; 
//   Rcpp::CharacterVector fea_name1 = adj_long["Var1"];
//   Rcpp::CharacterVector fea_name2 = adj_long["Var2"];
// 
//   List out(nrow_adj_long);
// 
//   for(int i = 0; i < nrow_adj_long - 1; i++){
// 
//     Rcpp::List df(3);
//     df.names() = Rcpp::CharacterVector::create("Feature_name1", "Feature_name2","NL_Names");
// 
//     CharacterVector feaNames1;
//     CharacterVector feaNames2;
//     CharacterVector NLNames;
// 
//     for(int j = 0; j < nrow_NL - 1; j++ )
//     {
//       float diff_mz1= std::fabs(NL_mz[j]-mz_diff[i]);
//      // float diff_mz2= std::fabs(Na_H_mz[j]-mz_diff[i]);
//      
//       if((diff_mz1<=diff_mz_thr || diff_mz1<=NL_mz[j]*ppm* 1e-6) && rt_diff[i]<=diff_rt_thr){
//            //||(diff_mz2<=diff_mz_thr||diff_mz2<=Na_H_mz[j]*ppm* 1e-6 && rt_diff[i]<=diff_rt_thr)
// 
//         feaNames1.push_back(fea_name1[i]);
//         feaNames2.push_back(fea_name2[i]);
//         NLNames.push_back(NL_name[j]);
//       }
//       
//     }
// 
//     df["Feature_name1"] = feaNames1;
//     df["Feature_name2"] = feaNames2;
//     df["NL_Names"]=NLNames;
//     out[i] = df;
// 
//   }
//   return out;
// 
// 
// }


#include <Rcpp.h>
using namespace Rcpp;

// NL match
// [[Rcpp::export]]
Rcpp::List findNL_fun(
    Rcpp::DataFrame NL_data,
    Rcpp::DataFrame adj_long,
    double diff_mz_thr,
    double diff_rt_thr,
    int ppm) {
  
  int nrow_adj_long = adj_long.nrow();
  int nrow_NL = NL_data.nrow();
   
  Rcpp::NumericVector mz_diff = adj_long["mz_diff"];
  Rcpp::NumericVector rt_diff = adj_long["rt_diff"];
  Rcpp::NumericVector NL_mz = NL_data["Accurate.Mass"];
  Rcpp::CharacterVector NL_name = NL_data["Neutral.Loss"];
  Rcpp::CharacterVector fea_name1 = adj_long["Var1"];
  Rcpp::CharacterVector fea_name2 = adj_long["Var2"];

  Rcpp::NumericVector mz1 = adj_long["mz_x"];  
  Rcpp::NumericVector mz2 = adj_long["mz_y"];  
  
  List out(nrow_adj_long);
   
  for (int i = 0; i < nrow_adj_long - 1; i++) {
    Rcpp::List df(3);
    df.names() = Rcpp::CharacterVector::create("Feature_name1", "Feature_name2", "NL_Names");
     
    CharacterVector feaNames1;
    CharacterVector feaNames2;
    CharacterVector NLNames;
     
    for (int j = 0; j < nrow_NL - 1; j++) {
      float diff_mz1 = std::fabs(NL_mz[j] - mz_diff[i]);
       
      if ((diff_mz1 <= diff_mz_thr || diff_mz1 <= NL_mz[j] * ppm * 1e-6) && rt_diff[i] <= diff_rt_thr) {
         
        
        if (mz1[i] < mz2[i]) {
          feaNames1.push_back(fea_name2[i]); // swap names
          feaNames2.push_back(fea_name1[i]);
         } else {
          feaNames1.push_back(fea_name1[i]);
          feaNames2.push_back(fea_name2[i]);
        }
         
        NLNames.push_back(NL_name[j]);
      }
    }
     
    df["Feature_name1"] = feaNames1;
    df["Feature_name2"] = feaNames2;
    df["NL_Names"] = NLNames;
    out[i] = df;
  }
   
  return out;
} 

