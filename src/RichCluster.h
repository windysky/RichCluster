#ifndef RICHCLUSTER_H
#define RICHCLUSTER_H

#include <Rcpp.h>
#include <string>
#include <functional>
#include "DistanceMetric.h"
#include "DistanceMatrix.h"
#include "AdjacencyList.h"
#include "LinkageMethod.h"
#include "ClusterList.h"

// class where the high level RichCluster algorithm is executed
class RichCluster {
public:
  static constexpr double SAME_TERM_DISTANCE = -99;

  RichCluster(Rcpp::CharacterVector termNameColumn,
              Rcpp::CharacterVector geneIDColumn,
              DistanceMetric dm, LinkageMethod lm);

  void computeDistances();
  void filterSeeds(); // informally denoting (node, neighbors) =: seed
  void mergeClusters();

  std::function<double(int, int)> distanceFunction() {
    return [this](int i, int j) {
      return _distanceMatrix.getDistance(i, j);
    };
  }

  // R export utilities
  Rcpp::NumericMatrix exportR_DistanceMatrix() {
    return distanceMatrix.export_RMatrix();
  }
  Rcpp::DataFrame exportR_ClusterList() {
    return clusterList.export_RDataFrame();
  }

private:
  std::unordered_set<int> filterSeed(int root, const std::unordered_set<int>& neighbors);
  ClusterList::ClusterIt RichCluster::findBestMergePartner(
      ClusterList::ClusterIt it1, std::list<std::unordered_set<int>>& clusters
  );

  // useful setup
  std::vector<std::string> _termNames;
  std::vector<std::string> _geneIDstrings; // each string contains all geneIDs for a single term
  int _nterms;
  // data structures
  DistanceMatrix distanceMatrix;
  AdjacencyList adjList;
  ClusterList clusterList;
  // distance/linkage classes
  DistanceMetric dm;
  LinkageMethod lm;
};

#endif
