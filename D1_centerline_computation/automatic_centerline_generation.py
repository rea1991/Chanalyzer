import numpy as np
import open3d as o3d
import os
import time
import vtk
from vtk.util.numpy_support import vtk_to_numpy
from vmtk import vmtkscripts
from vmtk import vtkvmtk

source  = [56.4254, 58.0546, 7.95225]
target  = [59.9422, 57.4528, 75.4727]

for C in range(1, 3):

    start = time.time()
    conformation = str(C).zfill(3)
    print(f"CONFORMATION {conformation}")
    if not os.path.exists("./output/frame." + conformation):
        os.mkdir("./output/frame." + conformation)

    mesh = o3d.io.read_triangle_mesh("../B_channel_projection/output/frame." + conformation + "/channel.off")
    mesh = o3d.geometry.TriangleMesh.compute_triangle_normals(mesh)
    o3d.io.write_triangle_mesh("../B_channel_projection/output/frame." + conformation + "/channel.stl", mesh)

    # Loading STL file:
    reader = vtk.vtkSTLReader()
    reader.SetFileName("../B_channel_projection/output/frame." + conformation + "/channel.stl")
    reader.Update()
    surface = reader.GetOutput()
    
    # Capping:
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = surface
    capper.Interactive = 0
    capper.Method = "centerpoint"
    capper.TriangleOutput = 0
    capper.Execute()
    surface = capper.Surface

    # Parameters for centerline computation:
    centerlines = vmtkscripts.vmtkCenterlines()
    centerlines.Surface = surface
    centerlines.Resampling = 1
    centerlines.ResamplingStepLength = 0.1
    centerlines.SeedSelectorName = "pointlist"
    centerlines.SourcePoints = source
    centerlines.TargetPoints = target
    
    # Centerline computation:
    centerlines.Execute()

    # Retrieving the centerline and some geometric properties:
    attributes = vmtkscripts.vmtkCenterlineAttributes()
    attributes.Centerlines = centerlines.Centerlines
    attributes.Execute()
    clNumpyAdaptor = vmtkscripts.vmtkCenterlinesToNumpy()
    clNumpyAdaptor.Centerlines = attributes.Centerlines
    clNumpyAdaptor.Execute()
    numpyCenterlines = clNumpyAdaptor.ArrayDict

    # Saving to txt:
    np.savetxt("output/frame." + conformation + "/points.txt", numpyCenterlines["Points"], delimiter=';')
    np.savetxt("output/frame." + conformation + "/maximumInscribedSphereRadius.txt", numpyCenterlines["PointData"]["MaximumInscribedSphereRadius"], delimiter=';')
    np.savetxt("output/frame." + conformation + "/abscissas.txt", numpyCenterlines["PointData"]["Abscissas"], delimiter=';')
    np.savetxt("output/frame." + conformation + "/parallelTransportNormals.txt", numpyCenterlines["PointData"]["ParallelTransportNormals"], delimiter=';')

    end = time.time()
    print(f"ELAPSED TIME: {end-start} seconds\n\n")
