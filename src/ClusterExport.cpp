#include <Rcpp.h>
#include "RichCluster.h"
#include "DistanceMetric.h"
#include "LinkageMethod.h"
#include <string>

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
Rcpp::List RichCluster(std::string distanceMetric, double distanceCutoff,
                       std::string linkageMethod, double linkageCutoff,
                       Rcpp::CharacterVector termNameColumn,
                       Rcpp::CharacterVector geneIDColumn) {

  RichCluster CM(termNameColumn,geneIDColumn);
  DistanceMetric DM(distanceMetric, distanceCutoff);

  CM.calculateDistanceScores(DM);

  /* MergeStrategy requires the following parameters:
      std::string mergeStrategy, (ex: "DAVID")
      double mergeCutoff, (0-1)
      std::string membershipStrategy, (ex: "DAVID")
      double membershipCutoff, (0-1)
  */
  LinkageMethod MS("DAVID", 0.5, "DAVID", membershipCutoff);
  CM.filterSeeds();
  Rcpp::DataFrame FilteredSeedMap = CM.exportR_SeedMap();
  CM.mergeSeeds(MS);
  // return Rcpp::List::create(
  //   Rcpp::_["DistanceMatrix"] = CM.exportR_DistanceMatrix(),
  //   Rcpp::_["SeedMap"] = CM.exportR_SeedMap(),
  //   Rcpp::_["FilteredSeeds"] = FilteredSeedMap,
  //   Rcpp::_["MergedSeeds"] = CM.exportR_ClusterList()
  // );
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
