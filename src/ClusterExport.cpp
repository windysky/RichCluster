// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "RichCluster.h"
#include "ClusterManager.h" // Added include
#include "DistanceMetric.h"
#include "LinkageMethod.h"
#include <string>

// [[Rcpp::export]]
Rcpp::List toyListExample() {
  // Create a simple list with two elements: x = 42, y = "hello"
  return Rcpp::List::create(
    Rcpp::_["x"] = 42,
    Rcpp::_["y"] = "hello"
  );
}

//' @name RichCluster
//'
//' This function clusters terms within an enrichment result with based on gene
//' similarity using specified distance metrics and merging strategies. It
//' integrates the following steps:
//' \enumerate{
//'   \item{Initializes a ClusterManager with input gene data.}
//'   \item{Calculates distance scores using the specified distance metric and cutoff.}
//'   \item{Filters seeds based on the provided merging strategy and membership cutoff.}
//'   \item{Merges seeds based on the merging strategy.}
//' }
//'
//' @param distanceMetric A string specifying the distance metric to use (e.g., "kappa").
//' @param distanceCutoff A double specifying the distance cutoff value.
//' @param mergeStrategy A string specifying the merge strategy to use (e.g., "DAVID").
//' @param membershipCutoff A double specifying the membership cutoff value (between 0 and 1).
//' @param termNameColumn A CharacterVector containing term names.
//' @param geneIDColumn A CharacterVector containing gene IDs.
//' @param PvalueColumn A NumericVector containing p-values.
//'
//' @return An R List object containing the following elements:
//' \itemize{
//'   \item{DistanceMatrix: The distance matrix used in clustering.}
//'   \item{SeedMap: The initial seed map of clusters.}
//'   \item{FilteredSeeds: The filtered seed map after applying the merge strategy.}
//'   \item{MergedSeeds: The final cluster list after merging seeds.}
//' }
//'
//' @examples
//' # Example usage
//' result <- RichCluster("kappa", 0.5, "DAVID", 0.7, termNames, geneIDs, pValues)
//' distanceMatrix <- result$distance_matrix
//' all_clusters <- result$all_clusters
//'
//' @export
// [[Rcpp::export]]
Rcpp::List RunRichCluster(std::string distanceMetric, double distanceCutoff,
                       std::string linkageMethod, double linkageCutoff,
                       Rcpp::CharacterVector termNameColumn,
                       Rcpp::CharacterVector geneIDColumn) {

  DistanceMetric dm_obj(distanceMetric, distanceCutoff); // Renamed DM to dm_obj for clarity
  LinkageMethod lm_obj(linkageMethod, linkageCutoff);   // Create LinkageMethod with correct constructor
  
  // Now “RichCluster” unambiguously refers to the class, not the function.
  RichCluster CM(termNameColumn, geneIDColumn, dm_obj, lm_obj);

  CM.computeDistances(); // Replaced calculateDistanceScores

  /* MergeStrategy requires the following parameters:
      std::string mergeStrategy, (ex: "DAVID")
      double mergeCutoff, (0-1)
      std::string membershipStrategy, (ex: "DAVID")
      double membershipCutoff, (0-1)
  */
  // LinkageMethod MS("DAVID", 0.5, "DAVID", membershipCutoff);
  CM.filterSeeds(); // Ensure this is present
  // Rcpp::DataFrame FilteredSeedMap = CM.exportR_SeedMap();
  CM.mergeClusters(); // Replaced mergeSeeds(MS)
  // return Rcpp::List::create(
  //   Rcpp::_["DistanceMatrix"] = CM.exportR_DistanceMatrix(),
  //   Rcpp::_["SeedMap"] = CM.exportR_SeedMap(),
  //   Rcpp::_["FilteredSeeds"] = FilteredSeedMap,
  //   Rcpp::_["MergedSeeds"] = CM.exportR_ClusterList()
  // );
  // Ensure the return statement matches the requirement (it already does from previous state)
  return Rcpp::List::create(
    Rcpp::_["distance_matrix"] = CM.exportR_DistanceMatrix(),
    Rcpp::_["all_clusters"] = CM.exportR_ClusterList()
  );
}

//' @name ComputeDistanceMatrix
//'
//' This function computes a distance matrix based on gene similarity using a specified distance metric.
//'
//' @param distanceMetric A string specifying the distance metric to use (e.g., "kappa").
//' @param distanceCutoff A double specifying the distance cutoff value.
//' @param termNameColumn A CharacterVector containing term names.
//' @param geneIDColumn A CharacterVector containing gene IDs.
//' @param PvalueColumn A NumericVector containing p-values.
//'
//' @return A NumericMatrix containing the distance matrix.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix ComputeDistanceMatrix(std::string distanceMetric, double distanceCutoff,
                                          Rcpp::CharacterVector termNameColumn,
                                          Rcpp::CharacterVector geneIDColumn,
                                          Rcpp::NumericVector PvalueColumn) {
  ClusterManager CM(termNameColumn, geneIDColumn, PvalueColumn);
  DistanceMetric DM(distanceMetric, distanceCutoff);

  CM.calculateDistanceScores(DM);
  return CM.exportR_DistanceMatrix();
}
