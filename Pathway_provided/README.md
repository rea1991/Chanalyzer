# Modified Pipeline: pathway is provided


### Before you start

It requires the following modules:
- itertools
- networkx
- numpy
- open3d
- os
- scipy
- sklearn
- sklearn
- time
- itertools
- math
- shutil


### Input files

- an off file representing the channel surface extracted by Step B of Chanalyzer;

- an xyzr file representing the pathway.


### Settings

- variables modname and selname (lines 168 and 169) has to be initialized as the input files (see above section);

- variable cutoff (line 170) is set by default to 3 and it represents the initial value of the maximum distance from the pathway for which a surface point of the channel is included into the reduced channel;

- variable mult (line 180) is set by default to 3 and its value is related to the halt condition for the while loop (line 182). Specifically, the algorithm will consider in the worst case mult*cutoff as the maximum distance from the pathway for which a surface point of the channel is included into the reduced channel.


### Output files

- an off file named 'channel_reduced.off' representing the reduced channel surface;

- a number of off files (the number depends on the while loop at line 182) named 'channel_reduced_i.off' representing some preliminary reduction of the channel surface ('channel_reduced_i.off' encodes the portion of the input channel surface at maximum distance i from the provided pathway).
