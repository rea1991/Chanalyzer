//
//  graph.h
//
//
//  Created by Andrea Raffo on 16/10/2020.
//

#ifndef graph_h
#define graph_h

#include <stdio.h>
#include <utility>
#include <vector>
#include <list>
#include <unordered_set>
#include <string>
#include <iostream>
#include <fstream>
#endif /* graph_h */

using namespace std;


#define infi 1000000000
using namespace std;

//CLASS FOR DIRECTED & WEIGHTED GRAPHS
class Graph {
private:
    // Number of vertices:
    int NumberOfVertices;
    
    // Properties (in the form of a vector of integers/doubles per vertex):
    std::vector<std::vector<int>>    properties_i;
    std::vector<std::vector<double>> properties_d;
    
    // Adjacency list (for each vertex, a list of pairs identifies children and weights):
    std::vector<std::list<pair<int, int>>> adj;
    
    // A recursive function used by DFS:
    void DFSUtil(int v, std::vector<bool>& visited,  std::vector<int>& flags, int K);
    
public:
    // Constructor:
    Graph(int NumberOfVertices);
    
    // Function to add an empty vertex to Graph:
    void addVertex();
    
    // Function to add an edge to Graph:
    void addEdge(int v, int w, int weight=1);
    
    // Functions to add a property (vector of integers or doubles) to a vertex v:
    void addPropertyI(int v, int prop);
    void addPropertyD(int v, double prop);
    
    // Get methods:
    int getNumberOfVertices();
    std::vector<int>    getPropertiesI(int v);
    std::vector<double> getPropertiesD(int v);
    std::vector<std::vector<int>>    getPropertiesI();
    std::vector<std::vector<double>> getPropertiesD();
    std::vector<std::list<pair<int, int>>>  getEdges();
    std::list<pair<int, int>>               getEdges(int v);
    
    //DIJKSTRA ALGORITHM FOR MINIMUM DISTANCE
    // Function to find the distance of the node from the given source vertex to the destination vertex:
    vector<int> dijkstraDist(int s, vector<int>& path);
    
    // DEPTH FIRST SEARCH (DFS) ALGORITHM
    // DFS traversal of the vertices reachable from v (only non-visited vertices are considered).
    // It uses recursive DFSUtil().
    void DFS(int v, std::vector<bool>& visited, std::vector<int>& flags, int K=0);
    
    //SAVING THE GRAPH TO TXT (UP TO NOW, JUST THE WEIGHTED GRAPH)
    void save_graph(std::string output_path);
    
    //LOADING THE GRAPH TO TXT (UP TO NOW, JUST THE WEIGHTED GRAPH)
    void load_graph(std::string input_path);

};
