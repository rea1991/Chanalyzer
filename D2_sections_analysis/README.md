# Step D2: sections' analysis


### Before you start

It requires the following modules:
- os
- math
- itertools
- numpy
- random
- networkx


### Input files

- an off file representing the channel surface extracted by the previous steps of Chanalyzer

- a txt file representing the sorted collection of vertices of the center line extracted by the previous steps of
  Chanalyzer (each line contains the 3 coordinates of a vertex)


### Settings

By varying the initialization of list ra (declared in line 247), one can choose the subset of center line points on which performing the section sections_analysis
Default: ra=range(normals.shape[0]); i.e., all the points of the center line are considered.


### Output files

- Center_Line_points.off; the input center line encoded as an off file

- Contour_Points_points.off; an off file representing the sections obtained by intersecting the planes orthogonal to the center line and the channel surface

- Contour_Points_points.txt; as above but encoded in a txt file

- Real_Contour_Points_points.off; an off file representing the VISIBLE PORTION of the sections obtained by intersecting the planes orthogonal to the center line and the channel surface

- Real_Contour_Points_points.txt; as above but encoded in a txt file

- Close_Points_points.off; an off file representing the collection of points P of the computed visible portions such that P is the closest to the center line

- Far_Points_points.off; an off file representing the collection of points P of the computed visible portions such that P is the farthest to the center line

- Min_Max_Distances_points.txt; a txt file representing, for each considered point P of the center line, the coordinates of P, the minimum distance between P and a point of the visible portion of the section corresponding to P, the maximum distance between P and a point of the visible portion of the section corresponding to P
