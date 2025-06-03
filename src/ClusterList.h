#ifndef CLUSTERLIST_H
#define CLUSTERLIST_H

#include "SeedMap.h"
#include "LinkageMethod.h"
#include "DistanceMatrix.h"
#include <Rcpp.h> // Added for Rcpp::DataFrame
#include <list>
#include <unordered_set>

class ClusterList{
public:
  using Cluster = std::unordered_set<int>;
  using ClusterIt = std::list<Cluster>::iterator;

  ClusterList() {};
  void addCluster(Cluster cluster) {clusterList.push_back(cluster);}
  void removeCluster(ClusterIt it) {clusterList.erase(it);}
  void mergeClusters(ClusterIt it1, ClusterIt it2) {
    if (it1 == it2) return;
    it1->insert(it2->begin(), it2->end());
    clusterList.erase(it2);
  }
  std::list<Cluster>& getList() {return clusterList;}
  Rcpp::DataFrame export_RDataFrame() { return Rcpp::DataFrame(); } // Added method

private:
  std::list<std::unordered_set<int>> clusterList;
};

#endif
