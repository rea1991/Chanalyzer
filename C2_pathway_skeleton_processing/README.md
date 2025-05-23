# Step C2_pathway: skeleton processing in case of a provided pathway


### Before you start

It requires the following modules:
- os
- math
- numpy
- networkx


### Input files

- The skeleton of a channel expressed as a graph extracted by Step C1 of Chanalyzer.
Each line represent an edge of the graph; e.g., line 2 55.1252 59.3976 76.5219 54.94 59.8402 76.9003 represents an edge whose vertices have coordinates (55.1252, 59.3976, 76.5219) and (54.94, 59.8402, 76.9003), respectively.


### Output files

- A file named Path_"name-of-input".off representing the path corresponding with the main direction in which the skeleton graph is developed. The output is expressed as a point cloud (subset of the nodes of the input graph) and it is encoded as an off file.

- As a secondary output, the program returns a file named Cloud_"name-of-input".off representing the nodes of the skeleton graph encoded as an off file.

The first and the last points inside Path_"name-of-input".off will be the source and target points in Step D1.
