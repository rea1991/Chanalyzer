import math
import numpy as np
import open3d as o3d
import os
import sys
import vtk
from vmtk import pypes
from vmtk import vmtkscripts
from vmtk import vtkvmtk

# Create a special environment to run VMTK
# Check http://www.vmtk.org/download/

# Parameters:
conformation = 1


# PREPROCESSING

# Output path:
conformation = str(conformation).zfill(3)
print(f"CONFORMATION {conformation}")
output_folder = "./output/frame." + conformation
if not os.path.exists(output_folder):
    os.mkdir(output_folder)
    
# Conversion to stl:
input_mesh_path = "../B_channel_projection/output/frame." + conformation
mesh = o3d.io.read_triangle_mesh(input_mesh_path + "/channel.off")
mesh = o3d.geometry.TriangleMesh.compute_triangle_normals(mesh)
o3d.io.write_triangle_mesh(input_mesh_path + "/channel.stl", mesh)



# VMTK AND GEOMETRIC INFORMATION

# Apply VMTK via the GUI:
command = "vmtkcenterlines -ifile " + input_mesh_path + "/channel.stl" + " -seedselector openprofiles -resampling 1 -resamplingstep 0.005 -ofile " + output_folder +  "/channel.vtp" + " --pipe vmtkcenterlineattributes -ifile " + output_folder +  "/channel.vtp" + " -ofile " + output_folder + "/channel.dat"
arg = pypes.PypeRun(command)

# Postprocess file:
points      = np.array([[0, 0, 0]])
MIRs        = np.array([0])
abscissas   = np.array([0])
PTNs        = np.array([[0, 0, 0]])

with open(output_folder + "/channel.dat") as lines:
    for n_line, line in enumerate(lines):
        if n_line>0:
            line_list   = line.split(" ")
            points      = np.vstack((points, np.array([[float(x) for x in line_list[0:3]]])))
            MIRs        = np.vstack((MIRs, np.array([float(line_list[3])])))
            abscissas   = np.vstack((abscissas, [float(line_list[4])]))
            PTNs        = np.vstack((PTNs, [[float(x) for x in line_list[5:8]]]))

np.savetxt(output_folder + "/points.txt", points[1:, :], delimiter=';')
np.savetxt(output_folder + "/MIRs.txt", MIRs[1:], delimiter=';')
np.savetxt(output_folder + "/abscissas.txt", abscissas[1:], delimiter=';')
np.savetxt(output_folder + "/PTNs.txt", PTNs[1:, :], delimiter=';')
