# Chanalyzer: a computational geometry approach for the analysis of protein channel shape and dynamics

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

### NEW RELEASE: main introduced features (please read the related PDF documentation)

Updates can be roughly classified into three main contributions:
- the introduction of a new step (Step E) aimed at providing new metrics for the retrieved channel and new visual representations of it;
- a modified pipeline (including new steps and programs) if a pathway is provided as input;
- a modified pipeline (including new steps and programs) to deal with occluded channels.

Please check the README files of the corresponding folders to learn more.

### Dataset behind the paper
The datasets analyzed for this study and the configuration parameter file for NanoShaper calculation can be found in the Zenodo repository, DOI: 10.5281/zenodo.6509652.

### Preprocessing tool
A script to convert a trajectory from the multiple pdb format to the frames that compose it in xyzr format is provided in https://github.com/concept-lab/mpdb2xyzr.git.

### References
[1]   A. Raffo, L. Gagliardi, U. Fugacci, L. Sagresti, S. Grandinetti, G. Brancato S. Biasotti, W. Rocchia. "Chanalyzer: a computational geometry approach for the analysis of protein channel shape and dynamics". *Frontiers in Molecular Biosciences* 9:933924, 2022. doi: 10.3389/fmolb.2022.933924

[2]   S. Decherchi, W. Rocchia. "A general and Robust Ray-Casting-Based Algorithm for Triangulating Surfaces at the Nanoscale". *PLOS ONE* 8(4): e59744, 2013.
