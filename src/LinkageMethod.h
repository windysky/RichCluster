#ifndef LINKAGEMETHOD_H
#define LINKAGEMETHOD_H

#include "DistanceMetric.h"
#include <Rcpp.h>
#include <unordered_set>
#include <string>
#include <functional>


class LinkageMethod {
public:
  LinkageMethod(std::string linkageMethod, double linkageThreshold)
    : _linkageMethod(linkageMethod), _linkageThreshold(linkageThreshold) {};

  double calculateLinkage(
      const std::unordered_set<int>& cluster1,
      const std::unordered_set<int>& cluster2
  );
  double getThreshold() {return _linkageThreshold;};
  void loadDistanceFunction(std::function<double(int, int)> distFct) {distanceFunction = distFct};

private:
  std::string _linkageMethod;
  double _linkageThreshold;
  std::function<double(int, int)> distanceFunction;

  // Types of supported membership strategies
  /*double david(const std::unordered_set<int>& cluster1,
              const std::unordered_set<int>& cluster2,
              std::function<double(int, int)> distanceFunction);*/
  double single(const std::unordered_set<int>& cluster1,
                const std::unordered_set<int>& cluster2);
  double complete(const std::unordered_set<int>& cluster1,
                  const std::unordered_set<int>& cluster2);
  double average(const std::unordered_set<int>& cluster1,
                const std::unordered_set<int>& cluster2);
  /*double ward(const std::unordered_set<int>& cluster1,
              const std::unordered_set<int>& cluster2,
              std::function<double(int, int)> distanceFunction);*/
};

#endif
