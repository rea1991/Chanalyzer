//
//  graph.cpp
//
//
//  Created by Andrea Raffo on 16/10/2020.
//

#include "Graph.h"


// Constructor:
Graph::Graph(int NumberOfVertices)
{
    this->NumberOfVertices = NumberOfVertices;
    adj = std::vector<std::list<pair<int, int>>>(NumberOfVertices);
    properties_i = std::vector<std::vector<int>>(NumberOfVertices);
    properties_d = std::vector<std::vector<double>>(NumberOfVertices);
}

// Function to add a vertex to Graph:
void Graph::addVertex(){
    NumberOfVertices++;
    adj.push_back({});
    properties_i.push_back({});
    properties_d.push_back({});
    
}

// Function to add an edge to Graph:
void Graph::addEdge(int v, int w, int weight)
{
    adj[v].push_back(std::pair<int, int>({w, weight})); // Add w to v’s list.
}

// Functions to add a property (vector of integers or doubles) to a vertex v:
void Graph::addPropertyI(int v, int prop)
{
    properties_i[v].push_back(prop);
}

void Graph::addPropertyD(int v, double prop)
{
    properties_d[v].push_back(prop);
}


// Get methods:
int Graph::getNumberOfVertices()
{
    return this->NumberOfVertices;
}

std::vector<int> Graph::getPropertiesI(int v)
{
    return this->properties_i[v];
}

std::vector<double> Graph::getPropertiesD(int v)
{
    return this->properties_d[v];
}

std::vector<std::vector<int>> Graph::getPropertiesI()
{
    return this->properties_i;
}
std::vector<std::vector<double>> Graph::getPropertiesD()
{
    return this->properties_d;
}

std::vector<std::list<pair<int, int>>> Graph::getEdges(){
    return this->adj;
}

std::list<pair<int, int>> Graph::getEdges(int v){
    return this->adj[v];
}


//DIJKSTRA ALGORITHM FOR MINIMUM DISTANCE
// Function to find the distance of the node from the given source vertex to the destination vertex:
vector<int> Graph::dijkstraDist(int s, vector<int>& path)
{
    // Stores distance of each vertex from source vertex:
    vector<int> dist(this->NumberOfVertices);

    // Boolean array that shows whether the vertex 'i' is visited or not:
    bool visited[this->NumberOfVertices];
    for (int i = 0; i < this->NumberOfVertices; i++) {
        visited[i] = false;
        path[i] = -1;
        dist[i] = infi;
    }
    dist[s] = 0;
    path[s] = -1;
    int current = s;

    // Set of vertices that has a parent (one or more) marked as visited:
    unordered_set<int> sett;
    while (true) {

        // Mark current as visited:
        visited[current] = true;
        std::list<pair<int, int>>::iterator i;
        for (i = adj[current].begin(); i != adj[current].end(); ++i){
            int v = (*i).first;
            if (visited[v])
                continue;

            // Inserting into the visited vertex:
            sett.insert(v);
            int alt = dist[current] + (*i).second;
            // Condition to check the distance is correct and update it if it is minimum from the previous computed distance
            if (alt < dist[v]) {
                dist[v] = alt;
                path[v] = current;
            }
        }
        sett.erase(current);
        if (sett.empty())
            break;

        // The new current
        int minDist = infi;
        int index = 0;

        // Loop to update the distance of the vertices of the Graph
        for (int a: sett) {
            if (dist[a] < minDist) {
                minDist = dist[a];
                index = a;
            }
        }
        current = index;
    }
    return dist;
}

// A recursive function used by DFS:
void Graph::DFSUtil(int v, std::vector<bool>& visited,  std::vector<int>& flags, int K)
{
    // Mark the current node as visited:
    visited[v] = true;
    flags[v] = K;
    
    // Recur for all the vertices adjacent to this vertex:
    std::list<pair<int, int>>::iterator i;
    for (i = adj[v].begin(); i != adj[v].end(); ++i)
    if (!visited[(*i).first])
        DFSUtil((*i).first, visited, flags, K);
}

// DFS traversal of the vertices reachable from v (only non-visited vertices are considered).
// It uses recursive DFSUtil().
void Graph::DFS(int v, std::vector<bool>& visited, std::vector<int>& flags, int K)
{
        // Call the recursive helper function to print DFS traversal:
        DFSUtil(v, visited, flags, K);
}


//SAVING THE GRAPH TO TXT (UP TO NOW, JUST THE WEIGHTED GRAPH)
void Graph::save_graph(std::string output_path){
    
    std::ofstream output;
    output.open(output_path);
    
    for (int current=0; current<adj.size();current++)
    {
        std::list<pair<int, int>>::iterator i;
        for (i = adj[current].begin(); i != adj[current].end(); ++i){
            int child_i = (*i).first;
            int weight_i = (*i).second;
            output << current << " " << child_i << " " << weight_i << std::endl;
        }
    }
    output.close();
}

//LOADING THE GRAPH TO TXT (UP TO NOW, JUST THE WEIGHTED GRAPH)
void Graph::load_graph(std::string input_path){
    std::ifstream input;
    input.open(input_path);
    int tmp1, tmp2, tmp3;
    while(input >> tmp1 >> tmp2 >> tmp3)
        addEdge(tmp1, tmp2, tmp3);
    input.close();
}
