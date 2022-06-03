# Chanalyzer: computational geometry for the analysis of protein channel shape and dynamics

This repository contains the routines behind the paper [1]. The repository is organized in four sub-folders:
- Step A:  extraction of the tetrahedral representation of the channel.
- Step B:  channel projection onto the SES generated via NanoShaper.
- Step C1: skeletonization.
- Step C2: extraction of source and target points.
- Step D1:  centerline computation and geometry processing.
- Step D2:  extraction and analysis of channel's sections.

A graphical abstract representing this pipeline is given in the following figure (note that C1 and C2 are summarized as a single step C; similarly, D1 and D2 are summed up as step D):

<p align="center">
<img src="https://github.com/rea1991/Chanalyzer/blob/a6468725c264985e9578e1e4075020e0588edba7/chanalyzer-core.png" alt="Chanalyzer-core" width="400"/>
</p>

### References
[1]   A. Raffo, L. Gagliardi, U. Fugacci, L. Sagresti, S. Grandinetti, G. Brancato S. Biasotti, W. Rocchia. "Chanalyzer: a computational geometry approach for the analysis of protein channel shape and dynamics", *Frontiers in Molecular Biosciences*, 2022.
