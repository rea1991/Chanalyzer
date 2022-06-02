import os
import math
import itertools
import numpy as np
import random
import networkx as nx


def LoadSurface(filename):

    infile = open("{}.off".format(filename), 'r')

    for line in itertools.islice(infile, 0, 2):
        pass

    nv, nt = int(line.split()[0]), int(line.split()[1])

    vertices = []
    edges = []
    triangles = []

    for i in range(nv):
        line = infile.readline()
        vertices.append([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])

    vertices=np.asarray(vertices)

    E = set()
    for line in itertools.islice(infile, 0, None):
        a = int(line.split()[1])
        b = int(line.split()[2])
        c = int(line.split()[3])

        edges = [[a, b], [a, c], [b, c]]
        for edge in edges:
            edge.sort()
            E.add(tuple(edge))


    edges = list(list(s) for s in E)

    edges=np.asarray(edges)

    print("\nLoaded Molecule {}".format(filename))
    print("# Vertices:{}".format(vertices.shape[0]))
    print("# Edges:{}".format(edges.shape[0]))

    infile.close()

    return vertices, edges

def LoadCenterLine(filename):

    path = os.path.realpath(__file__)
    dir = os.path.dirname(path)
    dir = dir.replace('D2_sections_analysis', 'D1_centerline_computation/output/frame.001/')

    infile = open(dir+"{}.txt".format(filename), 'r')

    points=[]

    for line in itertools.islice(infile, 0, None):

        new_point=[float(line.split(";")[0]),float(line.split(";")[1]),float(line.split(";")[2])]

        if points==[] or new_point != points[-1]:
             points.append(new_point)

    points=np.asarray(points)

    print("# Center Line Points:{}".format(points.shape[0]))

    f = open('output/Center_Line_{}.off'.format(filename), "w")
    f.write("OFF\n\n")
    f.write("{} 0 0 \n".format(points.shape[0]))

    for i in range(np.shape(points)[0]):
        f.write("{} {} {} \n".format(
            points[i, 0], points[i, 1], points[i, 2]))

    f.close()

    return points

def Int_Edge_Plane(p, n, p1, p2, v1, v2):

    D=p[0]*n[0]+p[1]*n[1]+p[2]*n[2]
    Int=[]

    num=n[0]*v1[0]+n[1]*v1[1]+n[2]*v1[2]-D
    den=n[0]*v1[0]-n[0]*v2[0]+n[1]*v1[1]-n[1]*v2[1]+n[2]*v1[2]-n[2]*v2[2]
    if den!=0:
        T=num/den
        if T>=0 and T<=1:
            q=list((1-T)*v1+T*v2)
            Int.append(q)

    return Int

def distance(q,p,r):


    den=(r[0]-p[0])*(r[0]-p[0])+(r[1]-p[1])*(r[1]-p[1])+(r[2]-p[2])*(r[2]-p[2])
    num=(r[0]-p[0])*q[0]+(r[1]-p[1])*q[1]+(r[2]-p[2])*q[2]-((r[0]-p[0])*p[0]+(r[1]-p[1])*p[1]+(r[2]-p[2])*p[2])

    t=num/den

    h=[p[0]+t*(r[0]-p[0]), p[1]+t*(r[1]-p[1]), p[2]+t*(r[2]-p[2])]

    cos=(r[0]-p[0])*(h[0]-p[0])+(r[1]-p[1])*(h[1]-p[1])+(r[2]-p[2])*(h[2]-p[2])

    if cos<=0:
        d=np.inf
    else:
        d=(q[0]-h[0])**2+(q[1]-h[1])**2+(q[2]-h[2])**2

    return d

def Closest_Point_WRT_Direction_R(intersection,p,r):

    points_with_distance=[]

    index=0
    for q in intersection:
        points_with_distance.append([q,distance(q,p,r),index])
        index=index+1

    points_with_distance=sorted(points_with_distance, key=lambda x: x[1])

    r_closets=intersection[0]
    r_index=0
    p_dist=np.inf
    l=max(int(len(points_with_distance)/35),10)

    for i in range(l):
        s=points_with_distance[i][0]
        d=(s[0]-p[0])**2+(s[1]-p[1])**2+(s[2]-p[2])**2
        if d<p_dist:
            r_closest=s
            r_index=points_with_distance[i][2]
            p_dist=d

    return r_closest, r_index

def Contour_Points(p,n,vertices,edges):

    num=n[0]*p[0]+n[1]*p[1]+n[2]*p[2]

    if n[0]!=0 and n[1]!=0:
        p1=[num/n[0], 0, 0]
        p2=[0, num/n[1], 0]

    elif n[0]!=0 and n[2]!=0:
        p1=[num/n[0], 0, 0]
        p2=[0, 0, num/n[2]]

    elif n[1]!=0 and n[2]!=0:
        p1=[0, num/n[1], 0]
        p2=[0, 0, num/n[2]]

    else:
        print("We have a problem!")

    intersection=[]

    for e in edges:
        Int=Int_Edge_Plane(p, n, p1, p2, vertices[e[0]], vertices[e[1]])

        for q in Int:
            intersection.append(q)

    adj_mat=np.zeros((len(intersection),len(intersection)))

    for a in range(len(intersection)):
        for b in range(a+1, len(intersection)):
            weight=math.sqrt((intersection[a][0]-intersection[b][0])**2+(intersection[a][1]-intersection[b][1])**2+(intersection[a][2]-intersection[b][2])**2)
            if weight<1:
                adj_mat[a][b]=1

    adj_mat=adj_mat+np.transpose(adj_mat)
    G = nx.from_numpy_matrix(adj_mat)

    H = nx.Graph()

    closest_points=[]

    for i in range(200):
        alpha=float(random.random()*np.pi)
        beta=np.arctan(-(n[2])/(n[0]*np.cos(alpha)+n[1]*np.sin(alpha)))
        rZero=np.sin(beta)*np.cos(alpha)+p[0]
        rOne=np.sin(beta)*np.sin(alpha)+p[1]
        rTwo=np.cos(beta)+p[2]
        r=[rZero,rOne,rTwo]
        r_closest, r_index = Closest_Point_WRT_Direction_R(intersection,p,r)
        closest_points.append(r_closest)
        H.add_nodes_from(nx.node_connected_component(G, r_index))

        beta=np.arctan(-(n[2])/(n[0]*np.cos(alpha)+n[1]*np.sin(alpha)))+np.pi
        rTwo=np.cos(beta)+p[2]
        rZero=np.sin(beta)*np.cos(alpha)+p[0]
        rOne=np.sin(beta)*np.sin(alpha)+p[1]
        r=[rZero,rOne,rTwo]
        r_closest, r_index =Closest_Point_WRT_Direction_R(intersection,p,r)
        closest_points.append(r_closest)
        H.add_nodes_from(nx.node_connected_component(G, r_index))

    closest_component=[]

    for u in H.nodes():
        closest_component.append(intersection[u])

    closest_points=closest_component

    min_dist=np.inf
    max_dist=0

    close=p
    far=p

    for w in closest_points:
        dist=math.sqrt((w[0]-p[0])**2+(w[1]-p[1])**2+(w[2]-p[2])**2)
        if dist<min_dist:
            close=w
            min_dist=dist
        if dist>max_dist:
            far=w
            max_dist=dist

    return intersection, closest_points, close, far, min_dist, max_dist


offname="pocket"
centerlinename="points"

vertices, edges = LoadSurface(offname)
points = LoadCenterLine(centerlinename)

normals = np.zeros((points.shape[0]-1, points.shape[1]))

for i in range(normals.shape[0]):
    normals[i,:]=points[i+1,:]-points[i,:]

cont_points=[]
real_cont_points=[]
close_points=[]
far_points=[]
close_dist=[]
far_dist=[]

ra=range(normals.shape[0])
# ra=range(211,212)

for i in ra:

    C,R,close,far,min_dist,max_dist=Contour_Points(points[i],normals[i],vertices,edges)
    cont_points.append(C)
    real_cont_points.append(R)
    close_points.append(close)
    far_points.append(far)
    close_dist.append(min_dist)
    far_dist.append(max_dist)
    print(i)

cont = sum( [ len(listElem) for listElem in cont_points])

real_cont = sum( [ len(listElem) for listElem in real_cont_points])

f = open('output/Contour_Points_{}.off'.format(centerlinename), "w")
f.write("OFF\n\n")
f.write("{} 0 0 \n".format(cont))

g = open('output/Contour_Points_{}.txt'.format(centerlinename), "w")
h = open('output/Min_Max_Distances_{}.txt'.format(centerlinename), "w")

f_close = open('output/Close_Points_{}.off'.format(centerlinename), "w")
f_close.write("OFF\n\n")
f_close.write("{} 0 0 \n".format(len(cont_points)))

f_far = open('output/Far_Points_{}.off'.format(centerlinename), "w")
f_far.write("OFF\n\n")
f_far.write("{} 0 0 \n".format(len(cont_points)))



for i in range(len(ra)):
    g.write("\n{} \n".format(ra[i]))
    g.write("{} {} {} \n".format(points[ra[i]][0],points[ra[i]][1],points[ra[i]][2]))
    h.write("{} {} {} {} {} \n".format(points[ra[i]][0],points[ra[i]][1],points[ra[i]][2],close_dist[i],far_dist[i]))
    f_close.write("{} {} {} \n".format(close_points[i][0],close_points[i][1],close_points[i][2]))
    f_far.write("{} {} {} \n".format(far_points[i][0],far_points[i][1],far_points[i][2]))

    for j in range(len(cont_points[i])):
        f.write("{} {} {} \n".format(cont_points[i][j][0], cont_points[i][j][1], cont_points[i][j][2]))
        g.write("{} {} {} \n".format(cont_points[i][j][0], cont_points[i][j][1], cont_points[i][j][2]))


f.close()
g.close()
h.close()
f_close.close()
f_far.close()

rf = open('output/Real_Contour_Points_{}.off'.format(centerlinename), "w")
rf.write("OFF\n\n")
rf.write("{} 0 0 \n".format(real_cont))

rg = open('output/Real_Contour_Points_{}.txt'.format(centerlinename), "w")

for i in range(len(ra)):
    rg.write("\n{} \n".format(ra[i]))
    rg.write("{} {} {} \n".format(points[ra[i]][0],points[ra[i]][1],points[ra[i]][2]))

    for j in range(len(real_cont_points[i])):
        rf.write("{} {} {} \n".format(real_cont_points[i][j][0], real_cont_points[i][j][1], real_cont_points[i][j][2]))
        rg.write("{} {} {} \n".format(real_cont_points[i][j][0], real_cont_points[i][j][1], real_cont_points[i][j][2]))

rf.close()
rg.close()
