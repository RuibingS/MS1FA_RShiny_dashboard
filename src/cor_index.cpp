#include <Rcpp.h>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
using namespace Rcpp;

// find the group index of a feature
std::string findGroup(std::unordered_map<std::string, std::string>& groups, 
                      const std::string& feature) {
  if (groups.find(feature) == groups.end() || groups[feature] == feature) {
    return feature;
  } else {
    return groups[feature] = findGroup(groups, groups[feature]);
  }
}

// check if merging two groups will satisfy the correlation threshold
bool MergeGroups(const std::unordered_map<std::string, std::unordered_map<std::string, double>>& correlationMap,
                 const std::unordered_set<std::string>& group1, const std::unordered_set<std::string>& group2, 
                 double threshold) {
  for (const std::string& f1 : group1) {
    for (const std::string& f2 : group2) {
      if (correlationMap.find(f1) == correlationMap.end() || 
          correlationMap.at(f1).find(f2) == correlationMap.at(f1).end() || 
          correlationMap.at(f1).at(f2) < threshold) {
        return false;
      }
    }
  } 
  return true;
} 

// [[Rcpp::export]]
List groupCorFeatures(DataFrame adj_long, double threshold) {
  std::vector<std::string> Var1 = as<std::vector<std::string>>(adj_long["Var1"]);
  std::vector<std::string> Var2 = as<std::vector<std::string>>(adj_long["Var2"]);
  Rcpp::NumericVector Freq = adj_long["Freq"];
  
  std::unordered_map<std::string, std::string> groups;
  std::unordered_map<std::string, std::unordered_set<std::string>> groupMembers;
  std::unordered_map<std::string, std::unordered_map<std::string, double>> correlationMap;
  
  // initialize correlation map 
  for (int i = 0; i < Var1.size(); ++i) {
    std::string feature1 = Var1[i];
    std::string feature2 = Var2[i];
    double correlation = Freq[i];
    
    correlationMap[feature1][feature2] = correlation;
    correlationMap[feature2][feature1] = correlation;
    
    if (groups.find(feature1) == groups.end()) {
      groups[feature1] = feature1;
      groupMembers[feature1].insert(feature1);
    } 
    if (groups.find(feature2) == groups.end()) {
      groups[feature2] = feature2;
      groupMembers[feature2].insert(feature2);
    }
  } 
  
  // each correlation pair
  for (int i = 0; i < Var1.size(); ++i) {
    std::string feature1 = Var1[i];
    std::string feature2 = Var2[i];
    double correlation = Freq[i];
    
    if (correlation >= threshold) {
      std::string group1 = findGroup(groups, feature1);
      std::string group2 = findGroup(groups, feature2);
      
      if (group1 != group2) {
        // Check if merging the groups satisfies the correlation threshold
        if (MergeGroups(correlationMap, groupMembers[group1], groupMembers[group2], threshold)) {
          // Merge the groups
          groups[group2] = group1;
          groupMembers[group1].insert(groupMembers[group2].begin(), groupMembers[group2].end());
          groupMembers.erase(group2);
        }
      }
    }
  } 
  
  // output with feature groups
  std::vector<std::string> output_features;
  std::vector<std::string> output_groups;
  std::unordered_map<std::string, std::string> groupNames;
  
  for (const auto& entry : groups) {
    std::string feature = entry.first;
    std::string group = findGroup(groups, feature);
    if (groupNames.find(group) == groupNames.end()) {
      groupNames[group] = "corgroup" + std::to_string(groupNames.size() + 1);
    } 
    
    output_features.push_back(feature);
    output_groups.push_back(groupNames[group]);
  } 
  
  DataFrame output = DataFrame::create(_["feature_name"] = output_features, _["cor_group"] = output_groups);
  
  // output with correlation pairs for each group
  std::unordered_map<std::string, std::vector<std::tuple<std::string, std::string, double>>> groupPairs;
  for (int i = 0; i < Var1.size(); ++i) {
    std::string feature1 = Var1[i];
    std::string feature2 = Var2[i];
    double correlation = Freq[i];
    
    std::string group1 = findGroup(groups, feature1);
    std::string group2 = findGroup(groups, feature2);
    
    if (group1 == group2) {
      groupPairs[groupNames[group1]].push_back(std::make_tuple(feature1, feature2, correlation));
    }
  } 
  
  // a list of data frames for each group
  List groupList;
  for (const auto& entry : groupPairs) {
    std::vector<std::string> pairVar1, pairVar2;
    std::vector<double> pairFreq;
    
    for (const auto& pair : entry.second) {
      pairVar1.push_back(std::get<0>(pair));
      pairVar2.push_back(std::get<1>(pair));
      pairFreq.push_back(std::get<2>(pair));
    } 
    
    DataFrame df = DataFrame::create(_["Var1"] = pairVar1, _["Var2"] = pairVar2, _["Freq"] = pairFreq);
    groupList[entry.first] = df;
  } 
  
  // return a list -  groups and the group detail list
  return List::create(_["feature_groups"] = output, _["group_details"] = groupList);
} 
