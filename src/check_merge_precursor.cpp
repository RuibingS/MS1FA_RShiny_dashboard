#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
// Function to merge precursor m/z values within a specified threshold
List mergePrecursorMZ(List spectraList1, List spectraList2, double mzThreshold, double rtThreshold) {
  List mergedMZList; // List to store merged precursor m/z values with names
  
  // Iterate over each spectrum in spectraList1
  for (int i = 0; i < spectraList1.size(); i++) {
    List spectrum1 = as<List>(spectraList1[i]); // Extract spectrum1
    
    double precursorMZ1 = as<double>(spectrum1["precursorMZ"]); // Get precursor m/z value of spectrum1
    double rtime1 = as<double>(spectrum1["rtime"]); // Get retention time of spectrum1
    
    bool matched = false; // Flag to indicate if a match is found
    
    // Iterate over each spectrum in spectraList2 to find a match for spectrum1
    for (int j = 0; j < spectraList2.size(); j++) {
      List spectrum2 = as<List>(spectraList2[j]); // Extract spectrum2
      
      double precursorMZ2 = as<double>(spectrum2["precursorMZ"]); // Get precursor m/z value of spectrum2
      double rtime2 = as<double>(spectrum2["rtime"]); // Get retention time of spectrum2
      
      // Check if the differences in precursor m/z and retention time are within the thresholds
      double mzDifference = std::abs(precursorMZ1 - precursorMZ2);
      double rtDifference = std::abs(rtime1 - rtime2);
      
      // If the differences are within the thresholds, merge the precursor m/z values
      if (mzDifference <= mzThreshold && rtDifference <= rtThreshold) {
        // Add the merged precursor m/z value to the list
        double mergedMZ = (precursorMZ1 + precursorMZ2) / 2; // Take the average as the merged value
        std::string name = std::to_string(i + 1); // Name for the merged m/z value
        mergedMZList[name] = mergedMZ; // Add the merged m/z value to the list with its name
        matched = true; // Set the flag to true to indicate a match is found
        break; // Exit the loop as the match is found
      }
    }
    
    // If no match is found, add precursorMZ1 to the list
    if (!matched) {
      std::string name = std::to_string(i + 1); // Name for precursorMZ1
      mergedMZList[name] = precursorMZ1; // Add precursorMZ1 to the list with its name
    }
  }
  
  return mergedMZList; // Return the list of merged precursor m/z values with names
}
