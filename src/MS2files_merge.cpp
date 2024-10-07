#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

List restructureData(DataFrame s4_object) {
  int n = s4_object.nrows();
  List output(n); // Output list to store restructured spectra
  
  // extract relevant columns from the S4 object
  NumericVector precursorMZ = s4_object["precursorMZ"];
  NumericVector rtime = s4_object["rtime"];
  NumericVector precursorIntensity = s4_object["precursorIntensity"];
  List mz_list = s4_object["mz"];
  List intensity_list = s4_object["intensity"];
  
  //
  for (int i = 0; i < n; ++i) {
    NumericVector mz_values = mz_list[i];
    NumericVector intensity_values = intensity_list[i];
    
    // Create a spectrum list with the extracted data
    List spectrum = List::create(Named("precursorMZ") = precursorMZ[i],
                                 Named("rtime") = rtime[i],
                                  Named("precursorIntensity") = precursorIntensity[i],
                                  Named("mz") = mz_values,
                                 Named("intensity") = intensity_values);
    
    output[i] = spectrum; // Add the spectrum to the output list
  }   
  
  return output; // Return the restructured data
}   


// calculate PPM difference between two m/z values
double calculatePPMDifference(double mz1, double mz2) {
  return std::abs(mz1 - mz2) / ((mz1 + mz2) / 2) * 1e6;
}
 

// concatenate two spectra while keeping the maximum intensity for each m/z value
List concatKeepMax(NumericVector mz1, NumericVector intensity1, 
                   NumericVector mz2, NumericVector intensity2) {
  // When both mz1 and mz2 are empty
  if (mz1.length() == 0 && mz2.length() == 0) {
    return List::create(Named("mz") = NumericVector(0), Named("intensity") = NumericVector(0));
  }
  
  // When mz1 and intensity1 are empty, return mz2 and intensity2
  if (mz1.length() == 0) {
    return List::create(Named("mz") = mz2, Named("intensity") = intensity2);
  }  
  
  // When mz2 and intensity2 are empty, return mz1 and intensity1
  if (mz2.length() == 0) {
    return List::create(Named("mz") = mz1, Named("intensity") = intensity1);
  }  
  
  std::vector<double> mz_combined;
  std::vector<double> intensity_combined;
   
  // Add all mz1 and intensity1 values to the combined vectors
  for (size_t i = 0; i < mz1.size(); i++) {
    mz_combined.push_back(mz1[i]);
    intensity_combined.push_back(intensity1[i]);
  }  
  
  // loop through mz2 to append or replace values in the combined vectors
  for (size_t j = 0; j < mz2.size(); j++) {
    bool found = false;
    for (size_t i = 0; i < mz_combined.size(); i++) {
      if (std::abs(mz_combined[i] - mz2[j]) <= 0.001) { // MZ values are considered the same
        found = true;
        // Update intensity if mz2's intensity is higher
        if (intensity2[j] > intensity_combined[i]) {
          intensity_combined[i] = intensity2[j];
        }  
        break;
      }
    }  
    // if mz2's value wasn't found in mz_combined, add it to the vectors
    if (!found) {
      mz_combined.push_back(mz2[j]);
      intensity_combined.push_back(intensity2[j]);
    }
  }
   
 
  NumericVector mz_combined_nv = wrap(mz_combined);
  NumericVector intensity_combined_nv = wrap(intensity_combined);
   
  return List::create(Named("mz") = mz_combined_nv, Named("intensity") = intensity_combined_nv);
}  




List mergeTwoSpectra(List spectraList1, 
                     List spectraList2, 
                     double mzThreshold,
                     double ppmThreshold, 
                     double rtThreshold) {
  List result;
  int jStart = 0;  // Initialize the starting index for the second list
  
  for (int i = 0; i < spectraList1.size(); i++) {
    List spectrum1 = as<List>(spectraList1[i]);
    double precursorMZ1 = as<double>(spectrum1["precursorMZ"]);
    double rtime1 = as<double>(spectrum1["rtime"]);
    double precursorIntensity1 = as<double>(spectrum1["precursorIntensity"]);
    
    NumericVector mz1 = spectrum1["mz"];
    NumericVector intensity1 = spectrum1["intensity"];
    bool matchFound = false;
    List bestMatch; // To hold the best matching spectrum
    
    for (int j = jStart; j < spectraList2.size(); j++) {
      List spectrum2 = as<List>(spectraList2[j]);
      double precursorMZ2 = as<double>(spectrum2["precursorMZ"]);
      double rtime2 = as<double>(spectrum2["rtime"]);
      double precursorIntensity2 = as<double>(spectrum2["precursorIntensity"]);
      
      if (std::abs(rtime1 - rtime2) > rtThreshold) {
        break; // No need to search further if the retention time difference exceeds the threshold
      } 
      
      double ppmDifference = calculatePPMDifference(precursorMZ1, precursorMZ2);
      if (ppmDifference <= ppmThreshold || std::abs(precursorMZ1 - precursorMZ2) <= mzThreshold) {
        // if precursor m/z values match within the threshold, merge the spectra
        NumericVector mz2 = spectrum2["mz"];
        NumericVector intensity2 = spectrum2["intensity"];
        List merged = concatKeepMax(mz1, intensity1, mz2, intensity2);
        
        // keep the spectrum with the higher precursor intensity
        if (precursorIntensity2 > precursorIntensity1) {
          bestMatch = List::create(Named("precursorMZ") = precursorMZ2,
                                   Named("rtime") = rtime2,
                                   Named("precursorIntensity") = precursorIntensity2,
                                   Named("mz") = merged["mz"],
                                                       Named("intensity") = merged["intensity"]);
        } else { 
          bestMatch = spectrum1; // Keep the original spectrum1
        } 
        
        matchFound = true;
        jStart = j + 1; 
        break; 
      } 
    }
    
    if (!matchFound) {
      result.push_back(spectrum1); // If no match found, keep the original spectrum1
    } else { 
      result.push_back(bestMatch); // Otherwise, add the merged spectrum to the result
    }
  } 
  
  return result;
} 



// [[Rcpp::export]]
List mergeMS2Data(List MS2_list, 
                  double mzThreshold,
                  double ppmThreshold,
                  double rtThreshold) {
  if (MS2_list.size() == 0) { // Return an empty list
   
    return List::create();
  }
  
 
  List merged = MS2_list[0];  
  
 
  for (int i = 1; i < MS2_list.size(); ++i) {
    List current = MS2_list[i];  
     
   
    merged = mergeTwoSpectra(merged, current, mzThreshold, ppmThreshold, rtThreshold);
  } 
   
  
  return merged;
}  



