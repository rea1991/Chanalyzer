from itertools import islice
import networkx as nx
import numpy as np
import open3d as o3d
import os
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors
from sklearn import preprocessing
import time
import itertools
import math
import shutil

def LoadSurface(modname):

    INfile = open("{}.off".format(modname), "r")

    for line in itertools.islice(INfile, 0, 2):
        pass


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

    print("\nLoaded Model {}".format(modname))
    print("# Vertices: {}".format(len(vertices)))
    print("# Edges: {}".format(len(edges)))
    print("# Triangles: {}".format(len(triangles)))

    INfile.close()

    return vertices, edges, triangles


def createReducedOFF(vertices,close_vertices,triangles,new_labels,modname,cutoff,done,original_modname):

    close_triangles=[]

    for t in triangles:
        a=int(new_labels[t[0]])
        b=int(new_labels[t[1]])
        c=int(new_labels[t[2]])
        if a*b*c>0.5:
            close_triangles.append([a-1, b-1, c-1])


    if done:
        OFFfile = open("{}_reduced.off".format(original_modname), "w")
    else:
        OFFfile = open("{}_reduced_{}.off".format(modname,cutoff), "w")

    OFFfile.write("OFF\n")
    OFFfile.write("\n")

    nv = len(close_vertices)
    nt = len(close_triangles)

    OFFfile.write("{} {} 0\n".format(nv,nt))

    for v in close_vertices:
        newline="{} {} {}\n".format(v[0], v[1], v[2])
        OFFfile.write(newline)

    for t in close_triangles:
        newline="3 {} {} {}\n".format(t[0], t[1], t[2])
        OFFfile.write(newline)


    OFFfile.close()

    return 0

def reduceB(modname,selname,cutoff):

    print("\nCreating reduced channel for cutoff={}...".format(cutoff))

    SELar = np.loadtxt("{}.xyzr".format(selname))

    vertices, edges, triangles = LoadSurface(modname)

    D = vertices+edges+triangles

    close_vertices=[]
    new_labels=np.zeros(len(vertices))

    for j in range(len(vertices)):
        v=vertices[j]
        for i in range(SELar.shape[0]):
            d=math.sqrt((v[0]-SELar[i,0])**2+(v[1]-SELar[i,1])**2+(v[2]-SELar[i,2])**2)
            if d<=cutoff:
                close_vertices.append(v)
                new_labels[j]=len(close_vertices)
                break

    createReducedOFF(vertices,close_vertices, triangles,new_labels,modname,cutoff,False,modname)

    return 0

def checkCutOff(modnameOne, modnameTwo, i_cutoff, modname):

    print("\nChecking cutoff={} by comparing Models {} and {}...".format(i_cutoff, modnameOne, modnameTwo))

    verticesOne, edgesOne, trianglesOne = LoadSurface(modnameOne)
    verticesTwo, edgesTwo, trianglesTwo = LoadSurface(modnameTwo)

    Gi = nx.Graph()
    Gi.add_edges_from(edgesTwo)

    V_three=[]
    for j in range(len(verticesOne)):
        for k in range(len(verticesTwo)):
            if verticesOne[j]==verticesTwo[k]:
                V_three.append(k)
                break

    CC=sorted(nx.connected_components(Gi), key=len, reverse=True)

    Ci=CC[0]

    # if len(CC[0])<5*len(CC[1]):
    #     Ci=set.union(CC[0],CC[1])

    out=True

    if len(set(V_three).difference(set(Ci)))>0:
        out=False

    else:
        print("\nWriting OFF file...")

        close_vertices=[]
        new_labels=np.zeros(len(verticesTwo))

        for v in Gi.nodes():
            if v in Ci:
                close_vertices.append(verticesTwo[v])
                new_labels[v]=len(close_vertices)


        createReducedOFF(verticesTwo,close_vertices,trianglesTwo,new_labels,modnameOne,i_cutoff,True,modname)

    return out

modname="channel"
selname="clc1_pathway_new_within3"
cutoff=3

start = time.time()

reduceB(modname,selname,cutoff)

modnameOne="{}_reduced_{}".format(modname,cutoff)
modnameTwo="{}_reduced_{}".format(modname,cutoff)

original_cutoff=cutoff
mult=3

while checkCutOff(modnameOne, modnameTwo, cutoff, modname)==False and cutoff<=mult*original_cutoff:
    cutoff=cutoff+1
    reduceB(modname,selname,cutoff)
    modnameTwo="{}_reduced_{}".format(modname,cutoff)

if cutoff>mult*original_cutoff:
    shutil.copy("{}_reduced_{}.off".format(modname, original_cutoff), "{}_reduced.off".format(modname))

end = time.time()
print(f"Elapsed time (in seconds): {end - start}")
