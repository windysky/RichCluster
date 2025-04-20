#include "MergeStrategy.h"

// currently, this computes david multiple linkage membership
// (separate into separate function in future)
double MergeStrategy::calculateMembership(
  const std::unordered_set<int>& clusterGroup,
  std::function<double(int, int)> distanceFunction
) {
  int totalPairs = 0;
  int significantClosePairs = 0;
  for (auto it1 = clusterGroup.begin(); it1 != clusterGroup.end(); ++it1)
  {
    for (auto it2 = std::next(it1); it2 != clusterGroup.end(); ++it2)
    {
      // dereference term iterators (we passed in reference to clusterGroup)
      int term1 = *it1;
      int term2 = *it2;
      // get distance (eg, kappa score) between terms
      double termDistance = distanceFunction(term1, term2);
      totalPairs++;
      if (termDistance >= _membershipCutoff) {
        significantClosePairs++;
      }
    }
  }
  if (totalPairs == 0) {
    // throw error
    return 0;
  }
  double membershipScore = significantClosePairs / totalPairs;
  return membershipScore;
}

/*
double MergeStrategy::singleLinkage(
  const std::unordered_set<int>& clusterGroup,
  std::function<double(int, int)> distanceFunction
) {
  double minDist = 100; // a default stupid value
  for (auto i = clusterGroup.begin(); i!= clusterGroup.end(); ++i){
    for (auto j = std::next(i); j!= clusterGroup.end(); ++j)
    { // dereference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distanceFunction(t1, t2);
      if (dist < minDist) {
        minDist = dist;
      }
    }
  }
  return minDist;
}

double MergeStrategy::completeMembership(
  const std::unordered_set<int>& clusterGroup,
  std::function<double(int, int)> distanceFunction
) {
  double maxDist = -100; // a default stupid value
  for (auto i = clusterGroup.begin(); i!= clusterGroup.end(); ++i){
    for (auto j = std::next(i); j!= clusterGroup.end(); ++j)
    { // dereference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distanceFunction(t1, t2);
      if (dist > maxDist) {
        maxDist = dist;
      }
    }
  }
  return maxDist;
}

double MergeStrategy::averageMembership(
  const std::unordered_set<int>& clusterGroup,
  std::function<double(int, int)> distanceFunction
) {
  int totalTermPairs = 0;
  double totalDist = 0;
  for (auto i = clusterGroup.begin(); i!= clusterGroup.end(); ++i){
    for (auto j = std::next(i); j!= clusterGroup.end(); ++j)
    { // dereference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distanceFunction(t1, t2);
      totalDist += dist;
      totalTermPairs++;
    }
  }
  return totalDist / totalTermPairs;
}
*/
