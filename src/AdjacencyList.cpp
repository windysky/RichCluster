#include "AdjacencyList.h"
#include "StringUtils.h"

// add a similar enough neighbor to a node's adjacency list
void AdjacencyList::addNeighbor(int node, int neighbor) {
  auto it = adjList.find(node);
  if (it != adjList.end()) { // if key term exists
    it->second.insert(neighbor);
  } else {
    adjList[node] = {neighbor}; // create new set for that key
  }
}

// returns true if neighbor exists in adjList[node]'s unordered_set
bool AdjacencyList::hasNeighbor(int node, int neighbor) {
  auto it = adjList.find(node);
  if (it == adjList.end()) {
    // node doesn't exist (bruh)
    return false;
  }
  const auto& neighbors = it->second;
  bool neighborExists = neighbors.find(neighbor) != neighbors.end(); // O(1) 😏
  return neighborExists;
}

