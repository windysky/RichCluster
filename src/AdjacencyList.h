#ifndef ADJACENCYLIST_H
#define ADJACENCYLIST_H

#include <Rcpp.h>
#include <vector>
#include <string>

class AdjacencyList {
public:
  AdjacencyList(int nterms) : _nterms(nterms) {}

  void addNeighbor(int node, int neighbor);
  bool hasNeighbor(int node, int neighbor);

  // return reference to adjList (used to initialize ClusterList)
  const std::unordered_map<int, std::unordered_set<int>>& getAdjList() const {
    return adjList;
  }

private:
  std::unordered_map<int, std::unordered_set<int>> adjList;
  int _nterms; // used to initialize size
};

#endif
