from itertools import islice
import networkx as nx
import numpy as np
import open3d as o3d
import os
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors
from sklearn import preprocessing
from sklearn.decomposition import PCA
import time
import itertools
import math
import shutil

#PARAMETERS

dir = "../Example - Input Pathway/"

structure="structure"
pathway="A"
selname="{}_pathway_{}".format(structure,pathway)

cutoff = 7

#########################

def LoadSurface(structure,dir):

    INfile = open(dir+"triangulatedSurf.off", "r")

    for line in itertools.islice(INfile, 0, 2):
        pass

    line = INfile.readline()

    if line[0].isnumeric()==False:
        line = INfile.readline()

    nv, nt = int(line.split()[0]), int(line.split()[1])

    vertices = []
    edges = []
    triangles = []

    for i in range(nv):
        line = INfile.readline()
        vertices.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])

    E = set()
    for i in range(nt):
        line = INfile.readline()
        a = int(line.split()[1])
        b = int(line.split()[2])
        c = int(line.split()[3])

        triangles.append([a, b, c])

        edges = [[a, b], [a, c], [b, c]]
        for edge in edges:
            edge.sort()
            E.add(tuple(edge))

    edges = list(list(s) for s in E)

    print("\nLoaded Model {}".format(structure))
    print("# Vertices: {}".format(len(vertices)))
    print("# Edges: {}".format(len(edges)))
    print("# Triangles: {}".format(len(triangles)))

    INfile.close()

    return vertices, edges, triangles

def createReducedOFF(vertices, cc_vertices, triangles, structure, cutoff, dir, pathway, num):
    """
    Generates an OFF file representing the sub-portion of the surface S
    with only the selected vertices from cc_vertices.
    """

    # Map original vertex indices to new indices in the reduced mesh
    vertex_map = {old_idx: new_idx for new_idx, old_idx in enumerate(cc_vertices)}

    # Filter triangles: keep only those whose vertices are all in cc_vertices
    reduced_triangles = []
    for t in triangles:
        if t[0] in vertex_map and t[1] in vertex_map and t[2] in vertex_map:
            reduced_triangles.append([vertex_map[t[0]], vertex_map[t[1]], vertex_map[t[2]]])

    # Open the OFF file for writing
    off_filename = "{}channel.off".format(dir)
    with open(off_filename, "w") as OFFfile:
        OFFfile.write("OFF\n\n")
        OFFfile.write("{} {} 0\n".format(len(cc_vertices), len(reduced_triangles)))

        # Write vertex positions
        for v_idx in cc_vertices:
            v = vertices[v_idx]
            OFFfile.write("{} {} {}\n".format(v[0], v[1], v[2]))

        # Write the triangle faces
        for t in reduced_triangles:
            OFFfile.write("3 {} {} {}\n".format(t[0], t[1], t[2]))

    print("Reduced OFF file saved as:", off_filename)
    return 0

def reduceB(structure,selname,cutoff,dir,pathway):

    print("\nCreating reduced channel for cutoff={}...".format(cutoff))

    SELar = np.loadtxt(dir+"{}.xyz".format(selname))

    f = open(dir+"{}.off".format(selname), "w")
    f.write("OFF\n\n")
    f.write("{} 0 0 \n".format(SELar.shape[0]))

    for i in range(SELar.shape[0]):
        f.write("{} {} {} \n".format(SELar[i, 0], SELar[i, 1], SELar[i, 2]))

    f.close()

    x = SELar[:,0]
    y = SELar[:,1]
    z = SELar[:,2]

    coords = np.array((x, y, z)).T

    pca = PCA(n_components=1)
    pca.fit(coords)
    direction_vector = pca.components_[0]
    origin = np.mean(coords, axis=0)

    NN=100

    g = open(dir+"{}_LINE.off".format(selname), "w")
    g.write("OFF\n\n")
    g.write("{} 0 0 \n".format(2*NN-1))

    for i in range(NN):
        t=0.5*i
        g.write("{} {} {} \n".format(origin[0]+direction_vector[0]*t,origin[1]+direction_vector[1]*t,origin[2]+direction_vector[2]*t))
        if t>0:
            t=-t
            g.write("{} {} {} \n".format(origin[0]+direction_vector[0]*t,origin[1]+direction_vector[1]*t,origin[2]+direction_vector[2]*t))

    g.close()

    h = open(dir+"{}_LINE_origin-direction.txt".format(selname), "w")

    h.write("{} {} {} \n".format(origin[0],origin[1],origin[2]))
    h.write("{} {} {} \n".format(direction_vector[0],direction_vector[1],direction_vector[2]))

    h.close()



    vertices, edges, triangles = LoadSurface(structure,dir)

    D = vertices+edges+triangles


    Gi = nx.Graph()
    Gi.add_edges_from(edges)

    for j in range(len(vertices)):
        # print(j)
        v=vertices[j]
        close=False
        for i in range(SELar.shape[0]):
            d=math.sqrt((v[0]-SELar[i,0])**2+(v[1]-SELar[i,1])**2+(v[2]-SELar[i,2])**2)
            if d<=cutoff:
                close=True
                break
        if close==False:
            Gi.remove_node(j)

    CC=sorted(nx.connected_components(Gi), key=len, reverse=True)


    createReducedOFF(vertices,list(CC[0]),triangles,structure,cutoff,dir,pathway,0)

    if len(CC[0])<5*len(CC[1]):
        createReducedOFF(vertices,list(CC[1]),triangles,structure,cutoff,dir,pathway,1)

    return 0

# start = time.time()

reduceB(structure,selname,cutoff,dir,pathway)

# end = time.time()
# print(f"Elapsed time (in seconds): {end - start}")
