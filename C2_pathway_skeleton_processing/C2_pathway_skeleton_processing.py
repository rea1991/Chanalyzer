import os
import math
import numpy as np
import networkx as nx

# PARAMETERS

dir = "../Example - Input Pathway/"

structure="structure"
pathway="A"
CC="0"
selname="{}_pathway_{}".format(structure,pathway)
filename = "CGAL_skel-poly.polylines_{}_CC{}".format(pathway,CC)
SELar = np.loadtxt(dir+"{}.xyz".format(selname))

#########################

def ListIndex(mylist, item):

    try:
        index = mylist.index(item)

    except ValueError:
        index = -1

    return index


def LoadGraph(filename, dir):

    infile = open(dir+"{}.txt".format(filename), 'r')

    data = np.loadtxt(infile)

    data = data[:, 1:]

    vertices = []

    nv = np.shape(data)[0]
    adMat = np.zeros((nv+1, nv+1))

    for i in range(nv):

        vOne = (data[i, 0], data[i, 1], data[i, 2])
        vTwo = (data[i, 3], data[i, 4], data[i, 5])

        indexOne = ListIndex(vertices, vOne)
        indexTwo = ListIndex(vertices, vTwo)

        if indexOne < 0:
            indexOne = len(vertices)
            vertices.append(vOne)

        if indexTwo < 0:
            indexTwo = len(vertices)
            vertices.append(vTwo)

        adMat[indexOne, indexTwo] = 1

    nv = len(vertices)

    adMat = adMat[:nv, :nv]

    adMat = adMat+adMat.transpose()

    G = nx.from_numpy_array(adMat)

    vertArray = np.zeros((len(vertices), 3))

    for i in range(len(vertices)):
        vertArray[i, 0] = vertices[i][0]
        vertArray[i, 1] = vertices[i][1]
        vertArray[i, 2] = vertices[i][2]

    return G, vertArray


def distance(u, v, w, vertices):

    l = vertices[w, 0]-vertices[v, 0]
    m = vertices[w, 1]-vertices[v, 1]
    n = vertices[w, 2]-vertices[v, 2]
    t = (l*(vertices[u, 0]-vertices[v, 0])+m*(vertices[u, 1]
         - vertices[v, 1])+n*(vertices[u, 2]-vertices[v, 2]))/(l**2+m**2+n**2)
    xH = vertices[v, 0]+t*l
    yH = vertices[v, 1]+t*m
    zH = vertices[v, 2]+t*n

    dist = math.sqrt((vertices[u, 0]-xH)**2
                     + (vertices[u, 1]-yH)**2+(vertices[u, 2]-zH)**2)

    return dist


def alignment(v, w, sp, vertices):

    alignment = 0

    i = 0
    for u in sp:
        if u != v and u != w:
            alignment = alignment+distance(u, v, w, vertices)
        i = i+1

    alignment = alignment/i

    return alignment

def pathway_distance(sp, SELar, vertices):

    path_dist = 0

    i = 0
    for u in sp:
        dist=np.infty
        for j in range(SELar.shape[0]):
            d=math.sqrt((vertices[u, 0]-SELar[j,0])**2+(vertices[u, 1]-SELar[j,1])**2+(vertices[u, 2]-SELar[j,2])**2)
            if d < dist:
                dist=d
        path_dist = path_dist+dist
        i = i+1

    path_dist = path_dist/i

    return path_dist


def score(v, w, sp, vertices, len_exponent, SELar):

    align = alignment(v, w, sp, vertices)

    path_dist = pathway_distance(sp, SELar, vertices)

    if align == 0:
        score = 0

    else:
        score = len(sp)/(align*path_dist**2)

    return score


G, vertices = LoadGraph(filename,dir)

print("Nodes:{}".format(G.number_of_nodes()))
print("Edges:{}".format(G.number_of_edges()))
print("\n")

f = open(dir + 'Cloud_{}_{}_CC{}.off'.format(structure,pathway,CC), "w")
f.write("OFF\n\n")
f.write("{} 0 0 \n".format(G.number_of_nodes()))

for i in range(np.shape(vertices)[0]):
    f.write("{} {} {} \n".format(
        vertices[i, 0], vertices[i, 1], vertices[i, 2]))

f.close()


min_len = nx.number_of_nodes(G)/4
len_exponent = 3
i = 0
j = 1
sp = nx.shortest_path(G, source=i, target=j)
max_score = score(i, j, sp, vertices, len_exponent, SELar)


all_paths = dict(nx.all_pairs_shortest_path(G))

print("All shortest paths computed")

for v in range(G.number_of_nodes()):
    print(v)
    for w in range(0, v):

        sp = all_paths[v][w]

        if len(sp) < min_len:
            scoreValue = 0
        else:
            scoreValue = score(v, w, sp, vertices, len_exponent, SELar)

        if scoreValue > max_score:
            max_score = scoreValue
            i = v
            j = w


sp = nx.shortest_path(G, source=i, target=j)

f = open(dir + 'Path_{}_{}_CC{}.off'.format(structure,pathway,CC), "w")
f.write("OFF\n\n")
f.write("{} 0 0 \n".format(len(sp)))

for i in sp:
    f.write("{} {} {} \n".format(
        vertices[i, 0], vertices[i, 1], vertices[i, 2]))

f.close()
