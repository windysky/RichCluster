#include "LinkageMethod.h"

// currently, this computes david multiple linkage membership
// (separate into separate function in future)
double LinkageMethod::calculateLinkage(
  const std::unordered_set<int>& cluster1,
  const std::unordered_set<int>& cluster2
) {
  // wrapper to call the actual function
  double linkageScore = 0;
  if (LinkageMethod::_linkageMethod == "david") { // (not actually david's implementation)
    // linkageScore = LinkageMethod::david(cluster1, cluster2, distanceFunction);
    linkageScore = LinkageMethod::single(cluster1, cluster2);
  } else if (LinkageMethod::_linkageMethod == "single") {
    linkageScore = LinkageMethod::single(cluster1, cluster2);
  } else if (LinkageMethod::_linkageMethod == "complete") {
    linkageScore = LinkageMethod::complete(cluster1, cluster2);
  } else if (LinkageMethod::_linkageMethod == "average") {
    linkageScore = LinkageMethod::average(cluster1, cluster2);
  } /* else if (LinkageMethod::_linkageMethod == "ward") {
    linkageScore = LinkageMethod::ward(cluster1, cluster2, distanceFunction);
  } */
  return linkageScore;
}


double LinkageMethod::single(
  const std::unordered_set<int>& cluster1,
  const std::unordered_set<int>& cluster2
) {
  double minDist = 100; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i){
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j)
    {
      if (i==j) {
        continue;
      }
      // de-reference the two term iterators
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

double LinkageMethod::complete(
    const std::unordered_set<int>& cluster1,
    const std::unordered_set<int>& cluster2
) {
  double maxDist = 0; // a default stupid value
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i){
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j)
    {
      if (i==j) {
        continue;
      }
      // de-reference the two term iterators
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

double LinkageMethod::average(
    const std::unordered_set<int>& cluster1,
    const std::unordered_set<int>& cluster2
) {
  double totalDist = 0;
  int n_terms = 0;
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i){
    for (auto j = cluster2.begin(); j!= cluster2.end(); ++j)
    {
      if (i==j) {
        continue;
      }
      // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      // get distance between them
      double dist = distanceFunction(t1, t2);
      totalDist += dist;
      n_terms ++;
    }
  }
  return totalDist/n_terms;
}


/*
double LinkageMethod::ward(
    const std::unordered_set<int>& cluster1,
    const std::unordered_set<int>& cluster2,
    std::function<double(int, int)> distanceFunction
) {

  double totalDist = 0; // a default stupid value
  int n_terms = 0;

  // INTRA-CLUSTER 1 DISTANCES
  for (auto i = cluster1.begin(); i!= cluster1.end(); ++i){
    for (auto j = cluster1.begin(); j!= cluster1.end(); ++j)
    { // de-reference the two term iterators
      int t1 = *i;
      int t2 = *j;
      if (t1==t2){ // ignore same term comparisons
        continue;
      }
      // get distance between them
      double dist = distanceFunction(t1, t2);
      totalDist += dist;
      n_terms ++;
    }
  }
  return totalDist/n_terms;
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
