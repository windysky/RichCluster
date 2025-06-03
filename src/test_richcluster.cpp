#include <Rcpp.h> // For Rcpp::CharacterVector
#include "RichCluster.h"    // Includes DistanceMetric.h, LinkageMethod.h, etc.
#include <vector>
#include <string>
#include <iostream> // For potential debug output, not strictly necessary for this test

int main() {
    Rcpp::Rcout << "Attempting to create mock objects..." << std::endl;

    // Create mock/default instances for RichCluster constructor
    Rcpp::CharacterVector mockTermNames;
    mockTermNames.push_back("TermA");
    mockTermNames.push_back("TermB");

    Rcpp::CharacterVector mockGeneIDs;
    mockGeneIDs.push_back("Gene1,Gene2");
    mockGeneIDs.push_back("Gene2,Gene3");

    Rcpp::Rcout << "Mock Rcpp::CharacterVectors created." << std::endl;

    // Assuming "kappa" and "average" are valid arguments based on previous context/files
    // and that DistanceMetric/LinkageMethod constructors don't do complex validation
    // that would require more specific mock data.
    DistanceMetric dm_obj("kappa", 0.5);
    Rcpp::Rcout << "DistanceMetric object created." << std::endl;

    LinkageMethod lm_obj("average", 0.5);
    Rcpp::Rcout << "LinkageMethod object created." << std::endl;

    Rcpp::Rcout << "Attempting to instantiate RichCluster..." << std::endl;
    try {
        RichCluster testCM(mockTermNames, mockGeneIDs, dm_obj, lm_obj);
        Rcpp::Rcout << "RichCluster object instantiated successfully." << std::endl;
    } catch (const std::exception& e) {
        Rcpp::Rcerr << "Exception caught during RichCluster instantiation: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        Rcpp::Rcerr << "Unknown exception caught during RichCluster instantiation." << std::endl;
        return 1;
    }

    return 0;
}
