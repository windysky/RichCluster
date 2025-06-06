// [[Rcpp::depends(Rcpp)]]
#include <Rcpp.h>
#include "RichCluster.h"
#include "ClusterManager.h" // Needed for ComputeDistanceMatrix
#include "DistanceMetric.h" // Needed for ComputeDistanceMatrix (and RichCluster via dm_obj)
#include "LinkageMethod.h"  // Needed for RichCluster (via lm_obj)
#include <string>

//' @name RichCluster
//' @export

// [[Rcpp::export]]
Rcpp::List RunRichCluster(std::string distanceMetric, double distanceCutoff,
                          std::string linkageMethod, double linkageCutoff,
                          Rcpp::CharacterVector termNameColumn,
                          Rcpp::CharacterVector geneIDColumn)  {

  // // Uncomment conversions from SEXP
  // std::string distanceMetric = Rcpp::as<std::string>(distanceMetricSEXP);
  // double distanceCutoff = Rcpp::as<double>(distanceCutoffSEXP);
  // std::string linkageMethod = Rcpp::as<std::string>(linkageMethodSEXP);
  // double linkageCutoff = Rcpp::as<double>(linkageCutoffSEXP);
  // Rcpp::CharacterVector termNameColumn = Rcpp::as<Rcpp::CharacterVector>(termNameColumnSEXP);
  // Rcpp::CharacterVector geneIDColumn = Rcpp::as<Rcpp::CharacterVector>(geneIDColumnSEXP);

  // Remove internal mocks
  // Rcpp::CharacterVector mockTermNames;
  // mockTermNames.push_back("MockTerm1");
  // Rcpp::CharacterVector mockGeneIDs;
  // mockGeneIDs.push_back("MockGene1");
  
  // Construct the helper objects exactly as the RichCluster constructor expects:
  DistanceMetric dm_obj(distanceMetric, distanceCutoff); // Use variables from SEXP parameters
  LinkageMethod lm_obj(linkageMethod, linkageCutoff); // Use variables from SEXP parameters

  // Now “RichCluster” refers unambiguously to the C++ class, not an R wrapper.
  RichCluster CM(termNameColumn, geneIDColumn, dm_obj, lm_obj); // Use variables from SEXP parameters

  CM.computeDistances();
  CM.filterSeeds();
  CM.mergeClusters();

  return Rcpp::List::create(
    Rcpp::_["distance_matrix"] = CM.exportR_DistanceMatrix(),
    Rcpp::_["all_clusters"] = CM.exportR_ClusterList()
  );
}

//' @name ComputeDistanceMatrix
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
