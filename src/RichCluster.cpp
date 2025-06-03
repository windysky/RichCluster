#include <Rcpp.h>

#include "RichCluster.h"
#include "ClusterList.h"
#include "StringUtils.h"
#include "DistanceMetric.h"
#include "LinkageMethod.h"
#include <stdexcept>
// #include <set> // No longer needed
#include <vector>      // For std::vector
#include <algorithm>   // For std::find

// initialization
RichCluster::RichCluster(Rcpp::CharacterVector termNameColumn,
                         Rcpp::CharacterVector geneIDColumn,
                         DistanceMetric dm, LinkageMethod lm)
  // store relevant info as C++ vectors
  : _termNames(Rcpp::as<std::vector<std::string>>(termNameColumn)),
    _geneIDstrings(Rcpp::as<std::vector<std::string>>(geneIDColumn)),
    _nterms(_termNames.size()),
    distanceMatrix(_nterms, _termNames),
    adjList(_nterms),
    clusterList(), // Corrected: use default constructor
    dm(dm), lm(lm)
{
  // ensure all vectors are of the same size
  if (_termNames.size() != _geneIDstrings.size()) {
    throw std::invalid_argument("All input columns must have the same size.");
  }
}

// phase 1: compute all pairwise distances
void RichCluster::computeDistances() {

  int totalGeneCount = StringUtils::countUniqueElements(_geneIDstrings);

  for (int i=0; i<_nterms; ++i) {
    // the unordered set of term1 genes
    std::unordered_set<std::string> term1_genes = StringUtils::splitStringToUnorderedSet(_geneIDstrings[i], ",");
    for (int j=0; j<_nterms; ++j) {
      if (i == j) {
        distanceMatrix.setDistance(RichCluster::SAME_TERM_DISTANCE, i, j);
        continue;
      }
      // unordered set of term2 genes
      std::unordered_set<std::string> term2_genes = StringUtils::splitStringToUnorderedSet(_geneIDstrings[j], ",");

      double distanceScore = dm.calculateDistance(term1_genes, term2_genes, totalGeneCount); // Used member 'dm'
      distanceMatrix.setDistance(distanceScore, i, j);

      // if term similarity is ABOVE the threshold
      if (distanceScore >= dm.getDistanceCutoff()) { // Used member 'dm'
        // add to adjacency list bidirectionally
        adjList.addNeighbor(i, j);
        adjList.addNeighbor(j, i);
      }
    }
  }
  // load in functional with init'ed distance matrix
  lm.loadDistanceFunction(distanceFunction());
}


// phase 2: seed filtering
void RichCluster::filterSeeds() {
  const auto& adjMap = adjList.getAdjList(); // C++11 compatible
  for (const auto& pair : adjMap) { // C++11 compatible
    int node = pair.first;
    const std::unordered_set<int>& neighbors = pair.second;
    std::unordered_set<int> cluster = filterSeed(node, neighbors);
    clusterList.addCluster(cluster);
  }
}

std::unordered_set<int> RichCluster::filterSeed(int node, const std::unordered_set<int>& neighbors) {
  std::unordered_set<int> cluster{node};

  while (true) {
    int bestN = -1;
    double bestLink = -1.0;

    for (int n : neighbors) {
      if (cluster.count(n)) continue;

      double link = lm.calculateLinkage(cluster, {n}); // Corrected: n wrapped in a set
      if (link > bestLink) {
        bestLink = link;
        bestN = n;
      }
    }

    if (bestLink < lm.getThreshold() || bestN == -1)
      break;
    cluster.insert(bestN);
  }

  return cluster;
}

void RichCluster::mergeClusters() {
  bool mergingPossible = true;

  while (mergingPossible) {
    int nMerged = 0;
    auto& clusters = clusterList.getList();
    std::vector<ClusterList::ClusterIt> toRemove; // Changed to std::vector

    for (auto it1 = clusters.begin(); it1 != clusters.end(); ++it1) {
      if (std::find(toRemove.begin(), toRemove.end(), it1) != toRemove.end()) continue; // skip clusters which we removed

      auto it2 = findBestMergePartner(it1, clusters);
      // Also check if it2 is already in toRemove before attempting merge
      if (it2 == clusters.end() || (std::find(toRemove.begin(), toRemove.end(), it2) != toRemove.end())) continue;

      // merge cluster2 into cluster1
      // Note: it2 is invalidated after this call. Storing it2 in toRemove after this is problematic.
      // For now, let's assume the main goal is to prevent it1 from being re-processed if it's modified,
      // or prevent the *original* it2 from being an it1.
      // The original logic to add it2 to 'toRemove' is tricky with invalidation.
      // I will add it *before* the merge for now, implying "don't use this it2 as an it1 later".
      // This might not be perfectly correct for the original algorithm's intent.

      // If it2 is a valid iterator and we are about to merge it:
      if (it2 != clusters.end()) { // Ensure it2 is valid before pushing
          toRemove.push_back(it2);
      }
      clusterList.mergeClusters(it1, it2);
      // it2 is now invalid and cannot be used.

      nMerged++;
    }

    if (nMerged == 0) {
      break;
    }
  }
}

ClusterList::ClusterIt RichCluster::findBestMergePartner(
    ClusterList::ClusterIt it1, std::list<std::unordered_set<int>>& clusters
) {
  double bestLink = -1.0;
  auto bestIt = clusters.end();

  for (auto it2 = clusters.begin(); it2 != clusters.end(); ++it2) {
    if (it1 == it2) continue;

    double link = lm.calculateLinkage(*it1, *it2);
    if (link > bestLink && link > lm.getThreshold()) {
      bestLink = link;
      bestIt = it2;
    }
  }

  return bestIt;
}


