#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
using namespace Rcpp;

// find the group index of a feature
std::string findGroup(std::unordered_map<std::string, std::string>& groups, 
                      const std::string& feature) {
  if (groups.find(feature) == groups.end() || groups[feature] == feature) {
    return feature;
  } else { 
    return groups[feature] = findGroup(groups, groups[feature]); // Path compression
  }
} 

// [[Rcpp::export]]
DataFrame groupFeatures(DataFrame df) {
  std::vector<std::string> ISF = as<std::vector<std::string>>(df["ISF"]);
  std::vector<std::string> PI = as<std::vector<std::string>>(df["PI"]);
  
  std::unordered_map<std::string, std::string> groups;
  std::unordered_map<std::string, std::string> groupNames;
  
  std::vector<std::string> uniqueFeatures(ISF.begin(), ISF.end());
  uniqueFeatures.insert(uniqueFeatures.end(), PI.begin(), PI.end());
  std::sort(uniqueFeatures.begin(), uniqueFeatures.end());
  uniqueFeatures.erase(std::unique(uniqueFeatures.begin(), uniqueFeatures.end()), uniqueFeatures.end());
  
  for (int i = 0; i < ISF.size(); i++) {
    std::string feature1 = ISF[i];
    std::string feature2 = PI[i];
    
    std::string group1 = findGroup(groups, feature1);
    std::string group2 = findGroup(groups, feature2);
    
    if (group1 != group2) {
      groups[group1] = group2;
    }
  }
  
  std::vector<std::string> output_features;
  std::vector<std::string> output_groups;
  
  for (const std::string& feature : uniqueFeatures) {
    std::string group = findGroup(groups, feature);
    
    if (groupNames.find(group) == groupNames.end()) {
      groupNames[group] =  "group" + std::to_string(groupNames.size() + 1);//
    }
    
    output_features.push_back(feature);
    output_groups.push_back(groupNames[group]);
  }
  
  DataFrame output = DataFrame::create(_["feature_name"] = output_features, _["group"] = output_groups);
  return output;
} 