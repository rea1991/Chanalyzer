import numpy as np
import math
import random
import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx

model = "Model"

threshold=1
sC=200
sC=1000

def FromCloudToOFF(vertices, offname):

    nv = vertices.shape[0]

    f = open('./output/{}.off'.format(offname), "w")
    f.write("OFF\n\n")
    f.write("{} 0 0 \n".format(nv))

    for i in range(nv):
        f.write("{} {} {} \n".format(
            vertices[i][0], vertices[i][1], vertices[i][2]))

    f.close()

    return 0


def createColoredOFFfunction(offname, vertices, Fv):

    nv=len(Fv)
    min_val = min(Fv)

    max_val=0

    for d in Fv:
         if ( d>max_val and d != np.inf):
            max_val=d

    check=0

    for d in Fv:
        if d == np.inf:
            d = max_val*2
            check=1

    if check == 1:
        max_val = max_val

    # use the coolwarm colormap that is built-in, and goes from blue to red
    cmap = mpl.cm.coolwarm
    normalized = mpl.colors.SymLogNorm(linthresh=1, vmin=min_val, vmax=max_val)

    # convert your distances to color coordinates
    Cv = cmap(normalized(Fv))

    f = open('./output/{}.off'.format(offname), "w")
    f.write("COFF\n\n")
    f.write("{} 0 0 \n".format(nv))

    for i in range(nv):
        f.write("{} {} {} {} {} {} {}\n".format(vertices[i][0], vertices[i][1], vertices[i][2], Cv[i,0],Cv[i,1],Cv[i,2],Cv[i,3]))

    f.close()

    return 0


def dist(u, v):

    return math.sqrt((u[0]-v[0])**2+(u[1]-v[1])**2+(u[2]-v[2])**2)

def closest_dist(u, U):
    Cdist=np.infty
    Idist=np.infty

    for i in range(np.shape(U)[0]):
        if dist(U[i], u) < Cdist:
            Cdist=dist(U[i], u)
            Idist=i

    return Cdist, Idist


def Length(data):

    l=0

    for i in range(np.shape(data)[0]-1):
        l=l+dist(data[i,:],data[i+1,:])

    return l

def distance(u, v, w):

    l = w[0]-v[0]
    m = w[1]-v[1]
    n = w[2]-v[2]
    t = (l*(u[0]-v[0])+m*(u[1]
         - v[1])+n*(u[2]-v[2]))/(l**2+m**2+n**2)
    xH = v[0]+t*l
    yH = v[1]+t*m
    zH = v[2]+t*n

    return dist(u, [xH,yH,zH])

def Turt(data):

    t = 0

    for i in range(1,np.shape(data)[0]-1):
        t = t+distance(data[0,:], data[i,:], data[np.shape(data)[0]-1,:])
        i = i+1

    t = t/np.shape(data)[0]

    return t

def Volume(data):

    vol=0

    hMax=0

    for i in range(data.shape[0]-1):

        h=dist(data[i], data[i+1])

        if h > hMax:
            hMax = h

        cone = (data[i][3]**2 + data[i][3]*data[i+1][3] + data[i+1][3]**2) * h

        vol = vol + cone

    vol = 1/3 * math.pi * vol

    return vol

### LOAD INPUT ###

infileC = "../D2_sections_analysis/output/Min_Max_Distances_points.txt"

chan = np.loadtxt(infileC)

nC = np.shape(chan)

print("Number of Vertices: {}".format(nC[0]))


### OFF FILES ###

createColoredOFFfunction(model+"_color", chan, chan[:,3])

### OFF FILES (SPHERES) ###

chan_spheres = []

for j in range(nC[0]):

    r=chan[j,3]

    for i in range(sC):
        A = random.uniform(0, 2*math.pi)
        z = random.uniform(-r,+r)

        v=[ math.sqrt(r**2-z**2)*math.cos(A)+chan[j,0], math.sqrt(r**2-z**2)*math.sin(A)+chan[j,1], z+chan[j,2], r ]

        toBeAddedC=True

        for k in range(nC[0]):
            if k!=j and dist(chan[k],v) <= r:
                toBeAddedC=False
                break
        if toBeAddedC==True:
            chan_spheres.append( v )

chan_spheres = np.asarray(chan_spheres, dtype=float)

createColoredOFFfunction(model+"_spheres", chan_spheres, chan_spheres[:,3])

# ### ANALYSIS ###

CXmap = [0]
CYmap = [chan[0,3]]

for i in range(1,nC[0]):

    CXmap.append(CXmap[i-1]+dist(chan[i-1,:],chan[i,:]))
    CYmap.append(chan[i,3])

plt.scatter(CXmap, CYmap, s=3, c='tab:blue')

plt.savefig("./output/Radius_{}.png".format(model))

print("\n- Length -")

print("Length: {}".format(round(Length(chan),2)))

print("\n- Straightness -")

print("Straightness: {}".format(round(1/Turt(chan),2)))

plt.show()
