# Step E: metrics and visualization


### Before you start

It requires the following modules:
- math
- numpy
- random
- networkx
- matplotlib


### Input files

- a txt file representing the centerline endowed with radius values of the maximal inscribed balls returned by Step D2.


### Settings

- Parameter `sC` sets the number of points adopted for sampling the spheres representing the channel returned by the previous steps of Chanalyzer.

- In order to run the program on models different than the one considered in the example modify line 144.


### Output files

Metrics about the centerline of the identified channel

- the number of vertices;
- the length;
- the straightness.

Informative visualizations of the channel:

- Model_color.off, the centerline of the retrieved channel colored in accordance with the radius values adopting the coolwarm colormap of Matplotlib (represented as a point cloud and encoded into an OFF file);

- Model_spheres.off, the retrieved channel (as a collection of spheres) colored in accordance with the radius values adopting the coolwarm colormap of Matplotlib (represented as a point cloud and encoded into an OFF file);

- Radius_Model.png, a PNG file representing the graph of the radius functions of the centerline of the considered model.



The tools and the code have been introduced in the repository

https://github.com/rea1991/GEO-Nav-methods

and discussed in [1].


### References

[1]   A. Raffo, U. Fugacci, S. Biasotti. "GEO-Nav: A geometric dataset of voltage-gated sodium channels", *Computers & Graphics*, vol. 115, pp. 285-295, 2023. DOI: 10.1016/j.cag.2023.06.023.
