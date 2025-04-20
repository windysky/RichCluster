#include "ClusterList.h"
#include "MergeStrategy.h"
#include "DistanceMatrix.h"
#include "StringUtils.h"
#include <functional>



void ClusterList::filterSeeds(MergeStrategy MS) {
  const auto& seedMapInstance = _seedMap.getSeedMap();

  for (const auto& termPair : seedMapInstance) {
    std::unordered_set<int> clusterGroup;
    clusterGroup.insert(termPair.first);
    clusterGroup.insert(termPair.second.begin(), termPair.second.end());

    double membership = MS.calculateMembership(clusterGroup, _distanceFunction);
    if (membership >= MS.getMembershipCutoff()) {
      clusterList.push_back(clusterGroup);
    }
  }
}


void ClusterList::mergeClusters(MergeStrategy MS) {
  bool mergingPossible = true;
  while (mergingPossible) {

    mergingPossible = false;
    int n_MergedGroups = 0; // stopping criteria, if 0 after outer loop then exit while loop

    // --- FIND BEST MERGE PARTNER LOOP ---
    auto it1 = clusterList.begin();
    while (it1 != clusterList.end()) {
      auto clusterGroup1 = *it1;

      int partnerIndex = 0;
      int bestMergePartnerIndex = -1; // index of mergePartner
      double bestMergeMembership = -99; // change to -inf later, (but mergeScore must be between 0-1 anyways)
      auto bestMergePartnerIt = clusterList.end();

      auto it2 = clusterList.begin();
      while (it2 != clusterList.end()) {
        if (it1 == it2) {
          ++it2;
          continue; // skip
        }

        // merge clusters and compute membership of merged group
        auto clusterGroup2 = *it2;
        std::unordered_set<int> mergedCluster;
        mergedCluster.insert(clusterGroup1.begin(), clusterGroup1.end());
        mergedCluster.insert(clusterGroup2.begin(), clusterGroup2.end());

        double mergeMembership = MS.calculateMembership(mergedCluster, _distanceFunction);

        // if above our desired cutoff
        if (mergeMembership >= MS.getMembershipCutoff()) {
          if (mergeMembership >= bestMergeMembership) {
            bestMergePartnerIndex = partnerIndex;
            bestMergeMembership = mergeMembership;
            bestMergePartnerIt = it2;
          }
        }
        ++it2;
        ++partnerIndex;
      }

      // --- THE MERGING, AFTER WE FOUND OUR PARTNER ---
      if (bestMergePartnerIndex != -1 && bestMergeMembership >= MS.getMembershipCutoff()) {
        std::unordered_set<int> mergedCluster;
        mergedCluster.insert(clusterGroup1.begin(), clusterGroup1.end());
        mergedCluster.insert(bestMergePartnerIt->begin(), bestMergePartnerIt->end());

        // remove original two clusters from clusterList
        it1 = clusterList.erase(it1);
        if (bestMergePartnerIt != clusterList.end()) {
          clusterList.erase(bestMergePartnerIt);
        }
        // and add the MERGED cluster back to clusterList
        clusterList.push_back(mergedCluster);

        n_MergedGroups++;
        mergingPossible = true;
        break; // Restart the loop after merging
      } else {
        ++it1;
      }
    }

    if (n_MergedGroups == 0) {
      return;
    }
  }
}

// export ClusterList as a R dataframe with termNames and index values
Rcpp::DataFrame ClusterList::export_RDataFrame() {
  std::vector<std::string> termIndicesColumn;
  std::vector<std::string> termNamesColumn;
  std::vector<int> clusterColumn;

  std::vector<std::string>& termNames = *(_seedMap._termNames); // dereference pointer

  int clusterNumber = 1;
  // Iterate through ClusterList and convert to vectors
  for (const auto& clusterGroup : clusterList) {
    // Turn the ints into a string
    std::string termIndicesString = StringUtils::unorderedSetToString(clusterGroup, ", ");
    termIndicesColumn.push_back(termIndicesString); // Append to termIndices

    std::vector<std::string> clusterGroupTerms;
    for (const auto& term_index : clusterGroup) {
      std::string term = termNames[term_index];
      clusterGroupTerms.push_back(term);
    }
    // Convert vector/unordered_set to one comma-delimited string
    std::string clusterGroupTerms_string = StringUtils::vectorToString(clusterGroupTerms, ", ");
    termNamesColumn.push_back(clusterGroupTerms_string); // Append to termNames

    // Append cluster number to clusterColumn
    clusterColumn.push_back(clusterNumber++);
  }

  // Create and return a DataFrame using Rcpp
  return Rcpp::DataFrame::create(Rcpp::Named("Cluster") = clusterColumn,
                                 Rcpp::Named("TermNames") = termNamesColumn,
                                 Rcpp::Named("TermIndices") = termIndicesColumn);
}
