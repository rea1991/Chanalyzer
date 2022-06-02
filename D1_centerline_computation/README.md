# Step D1: centerline extraction


### Before you start
This python implementation is tasked to load the channel extracted in Step B (in OFF format). The following external modules are required: `vtk` and `vmtk`.

### Input files
It requires a triangulated channel (as generated in Step B). In its automatic implementation, it requires a source and a target point, which is what is obtained in Step C2. If, for numerical instability reasons, the centerline is not satisfactory, then the user can either:
- Try to move slightly source and target points.
- Recompute source and target points by Step C2, whether these points were inherited from a previous conformation (and the conformations were aligned).  
- Use the semiautomatic implementation: it calls the VMTK GUI and let the user pick source and target points.

### Output files
It returns four TXT files:
- The centerline as a collection of points, one per line.
- The maximum inscribed sphere radii, one per point of the centerline.
- The parallel transport normals, one per point of the centerline.
- The curvilinear abscissas, one per point of the centerline.
