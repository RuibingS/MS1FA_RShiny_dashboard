#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
CharacterVector Check_Iso_Charge(CharacterVector name,
                                 NumericVector mz,
                                 NumericVector rt,
                                 NumericVector intensity,
                                 double iso_threshold = 1.003355, //  C13 isotope difference
                                 int max_charge = 3,
                                 int max_isotope = 4) {
  
  CharacterVector labels(name.size()); 
  int groupIndex = 1;
  
  for (int i = 0; i < name.size(); i++) {
    if (labels[i] != "") continue;
    
    std::vector<int> maxIsotopeForCharge(max_charge + 1, 0); 
    
    for (int j = i + 1; j < name.size(); j++) {
      if (labels[j] != "") continue;
      
      double mz_diff = mz[j] - mz[i];
      double rt_diff = std::abs(rt[j] - rt[i]);
      

      for (int c = 1; c <= max_charge; c++) {
        for (int m = 1; m <= max_isotope; m++) {
          double expected_mz_diff = m * iso_threshold / c;  
          double max_int_ratio = (mz[i] * c / 12) * 0.011;  // intensity ratio check
          
          if (m == 1) { // the first isotope [M+1]
            if (std::abs(mz_diff - expected_mz_diff) <= 0.002 && rt_diff < 3 && intensity[i] * max_int_ratio >= intensity[j]) {
            
              labels[i] = Rcpp::String("[" + std::to_string(groupIndex) + "][M]" +
                (c > 1 ? std::to_string(c) : "") + "+");
              
              labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                std::to_string(m) + "]" + (c > 1 ? std::to_string(c) : "") + "+");
              
              maxIsotopeForCharge[c] = m;  
            } 
          } else {  // for  isotopes M+2, M+3.. S34 will be annotated as M+2
            // intensity check to the monoisotopic peak 
            if (std::abs(mz_diff - expected_mz_diff) <= 0.01 && rt_diff < 3 && intensity[j] <= intensity[i]) {
              if (maxIsotopeForCharge[c] == m - 1) {  
                labels[j] = Rcpp::String("[" + std::to_string(groupIndex) + "][M+" +
                  std::to_string(m) + "]" + (c > 1 ? std::to_string(c) : "") + "+");
                
                maxIsotopeForCharge[c] = m;  
              }
            }
          }
        }
      }
    } 
    
    if (!labels[i].empty()) {
      groupIndex++;  
    } 
  }
  
  return labels;
} 
