from itertools import islice
import networkx as nx
import numpy as np
import open3d as o3d
import os
from scipy.sparse import coo_matrix
from sklearn.neighbors import NearestNeighbors
from sklearn import preprocessing
import time



# PARAMETERS
K = 20      # k for k-nearest neighbour
conformation = 2
root            = "../A_tetrahedral_representation/example/IonicChannels/mscl/5ligand"


conformation = str(conformation).zfill(3)
print(f"CONFORMATION {conformation}")
root = root + "/frame." + conformation + "/"
input_folder    = root
output_folder   = "./output/frame." + conformation + "/"

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


# FUNCTION TO PROCESS CONNECTED COMPONENTS
def tri_adiac_cell(mat_ver, mat_tri):

    num_ver = mat_ver.shape[0]
    num_tri = mat_tri.shape[0]
    
    succ = [1, 2, 0]
    t = [i for i in range(num_tri)]
    t_t = coo_matrix((num_ver, num_ver), dtype=np.int8)
    for i in range(3):
        t_t = t_t + coo_matrix((t, (mat_tri[:, i], mat_tri[:, succ[i]])), shape=(num_ver, num_ver))
    
    t_t = t_t.tocsr()
    mat_tri_bu = mat_tri.copy()
    for i in range(3):
        mat_tri = np.column_stack((mat_tri, t_t[mat_tri_bu[:, succ[i]], mat_tri_bu[:, i]].transpose()))
    
    ext_mat_tri = mat_tri.copy()
    mat_tri = mat_tri[:, 3:6]
    
    ad_tri = np.zeros(mat_tri.shape)
    
    for i in range(mat_tri.shape[0]):
        bb = []
        cc = 0
        for ii in range(3):
            if mat_tri[i, ii]!=-1:
                cc = cc+1
                bb = bb + [mat_tri[i, ii]]
        cc = 0
        ad_tri[i, :] = bb
    
    return (ad_tri.astype(int),ext_mat_tri)



# DATA LOADING AND PREPROCESSING

start = time.time()

# Load tetrahedral connected component:
T = np.loadtxt(root + "ED/ED_tet/CC_0.txt")

# Load SES triangle mesh:
triangulatedSurf = open(input_folder + "triangulatedSurf.off")
lines=triangulatedSurf.readlines()
triangulatedSurf.close()
num_of_simplices = [int(x) for x in lines[3].split()]

# Extract SES vertices and facets:
V = np.array([[float(x) for x in lines[i].strip().split(' ')] for i in range(4, 4 + num_of_simplices[0])])
F = np.array([[int(x) for x in lines[i].strip().split(' ')] for i in range(4 + num_of_simplices[0], 4 + num_of_simplices[0] + num_of_simplices[1])])

# Barycenters of tetrahedra (of the CCs) and of facets (of the SES)
B_tet_CC  = (T[:,0:3] + T[:,3:6] + T[:,6:9] + T[:,9:12])/4;
B_fac_SES = (V[F[:, 1], :] + V[F[:, 2], :] + V[F[:, 3], :])/3;

# k-nearest neighbours:
knn = NearestNeighbors(n_neighbors=K).fit(B_tet_CC)
distances, indices = knn.kneighbors(B_fac_SES)

# First facet of tetrahedra:
F1 = T[:, [0,1,2,3,4,5,9,10,11]]
# Second facet of tetrahedra:
F2 = T[:, [3,4,5,6,7,8,9,10,11]]
# Third facet of tetrahedra:
F3 = T[:, [0,1,2,6,7,8,3,4,5]]
# Fourth facet of tetrahedra:
F4 = T[:, [0,1,2,9,10,11,6,7,8]]

# Normals to the four facets of tetrahedra:
N_facets_F1 = preprocessing.normalize(np.cross(F1[:,3:6] - F1[:,0:3], F1[:,6:9] - F1[:,0:3]))
N_facets_F2 = preprocessing.normalize(np.cross(F2[:,3:6] - F2[:,0:3], F2[:,6:9] - F2[:,0:3]))
N_facets_F3 = preprocessing.normalize(np.cross(F3[:,3:6] - F3[:,0:3], F3[:,6:9] - F3[:,0:3]))
N_facets_F4 = preprocessing.normalize(np.cross(F4[:,3:6] - F4[:,0:3], F4[:,6:9] - F4[:,0:3]))

# Baricenters of the facets:
B_F1 = np.column_stack((
    (F1[:, 0] + F1[:, 3] + F1[:, 6])/3,
    (F1[:, 1] + F1[:, 4] + F1[:, 7])/3,
    (F1[:, 2] + F1[:, 5] + F1[:, 8])/3))
B_F2 = np.column_stack((
    (F2[:, 0] + F2[:, 3] + F2[:, 6])/3,
    (F2[:, 1] + F2[:, 4] + F2[:, 7])/3,
    (F2[:, 2] + F2[:, 5] + F2[:, 8])/3))
B_F3 = np.column_stack((
    (F3[:, 0] + F3[:, 3] + F3[:, 6])/3,
    (F3[:, 1] + F3[:, 4] + F3[:, 7])/3,
    (F3[:, 2] + F3[:, 5] + F3[:, 8])/3))
B_F4 = np.column_stack((
    (F4[:, 0] + F4[:, 3] + F4[:, 6])/3,
    (F4[:, 1] + F4[:, 4] + F4[:, 7])/3,
    (F4[:, 2] + F4[:, 5] + F4[:, 8])/3))

# Initializing 1d-array for colours (1 means that it must be coloured):
colors = np.zeros((B_fac_SES.shape[0], 1))

# Check which of the SES facets are "inside" a tetrahedron and color them:
for l in range(K):
    
    c1 = np.einsum("ij,ij->i", B_fac_SES - B_F1[indices[:, l], :], N_facets_F1[indices[:, l], :])
    c2 = np.einsum("ij,ij->i", B_fac_SES - B_F2[indices[:, l], :], N_facets_F2[indices[:, l], :])
    c3 = np.einsum("ij,ij->i", B_fac_SES - B_F3[indices[:, l], :], N_facets_F3[indices[:, l], :])
    c4 = np.einsum("ij,ij->i", B_fac_SES - B_F4[indices[:, l], :], N_facets_F4[indices[:, l], :])

    mask = np.all([c1<0, c2<0, c3<0, c4<0], axis=0)
    colors[np.where(mask)[0]] = 1



# CONNECTIVITY OF THE SES MESH

ad_tri, ext_mat_tri = tri_adiac_cell(V, F[:,1:4])
col_ad_tri = np.squeeze(colors[ad_tri])

# First edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 0][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 0][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges = np.hstack((a1[:, np.newaxis], ad_tri[a1, 0][:, np.newaxis]))
edges = np.vstack((edges,
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 0][:, np.newaxis]))
    ))

# Second edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 1][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 1][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges = np.vstack((edges,
    np.hstack((a1[:, np.newaxis], ad_tri[a1, 1][:, np.newaxis])),
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 1][:, np.newaxis]))
    ))
    
# Third edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 2][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri[:, 2][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges = np.vstack((edges,
    np.hstack((a1[:, np.newaxis], ad_tri[a1, 2][:, np.newaxis])),
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 2][:, np.newaxis]))
    ))
    
# Recover isolated vertices:
edges = np.vstack((edges, np.hstack((
    np.arange(0, F.shape[0], 1, dtype=int)[:, np.newaxis],
    np.arange(0, F.shape[0], 1, dtype=int)[:, np.newaxis]))))

# Define connectivity graph:
G = nx.from_edgelist(edges)
CCs         = [G.subgraph(c) for c in nx.connected_components(G)]
size_CCs    = [int(CCs[i].number_of_nodes()) for i in range(len(CCs))]



# REMOVE POSSIBLE IMPERFECTIONS

# Color only the biggest connected component:
sort_index = np.argsort(size_CCs)
for i in range(len(CCs)):
    if i!=sort_index[-1]:
        colors[CCs[i].nodes()] = 1

facets_gray = np.where(colors == 0)[0]
facets_red  = np.where(colors == 1)[0]
F_colored = np.hstack((F, np.ones((F.shape[0], 4))))
F_colored[facets_gray, 4:7] = 200*F_colored[facets_gray, 4:7]
F_colored[facets_red, 4]    = 255*F_colored[facets_red, 4]
F_colored[facets_red, 5]    = 0
F_colored[facets_red, 6]    = 0
F_colored[:, 7] = 255
col_ad_tri_new = np.squeeze(colors[ad_tri])

# First edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 0][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 0][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges_new = np.hstack((a1[:, np.newaxis], ad_tri[a1, 0][:, np.newaxis]))
edges_new = np.vstack((edges_new,
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 0][:, np.newaxis]))
    ))

# Second edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 1][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 1][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges_new = np.vstack((edges_new,
    np.hstack((a1[:, np.newaxis], ad_tri[a1, 1][:, np.newaxis])),
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 1][:, np.newaxis]))
    ))
    
# Third edge of SES facets:
a1 = np.where(np.all([
     (colors == False).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 2][:, np.newaxis] == False).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
a2 = np.where(np.all([
    (colors == True).any(axis=1)[:, np.newaxis],
     (col_ad_tri_new[:, 2][:, np.newaxis] == True).any(axis=1)[:, np.newaxis]
     ], axis=0))[0]
edges_new = np.vstack((edges_new,
    np.hstack((a1[:, np.newaxis], ad_tri[a1, 2][:, np.newaxis])),
    np.hstack((a2[:, np.newaxis], ad_tri[a2, 2][:, np.newaxis]))
    ))

# Recover isolated vertices:
edges_new = np.vstack((edges_new, np.hstack((
    np.arange(0, F.shape[0], 1, dtype=int)[:, np.newaxis],
    np.arange(0, F.shape[0], 1, dtype=int)[:, np.newaxis]))))

# New graph:
G_new = nx.from_edgelist(edges_new)
CCs_new      = [G_new.subgraph(c) for c in nx.connected_components(G_new)]
size_CCs_new = [int(CCs_new[i].number_of_nodes()) for i in range(len(CCs_new))]

# Filling possible holes:
sort_index_new = np.argsort(size_CCs_new)

for i in range(len(CCs_new)):
    if i!=sort_index_new[-1] and i!=sort_index_new[-2]:
        colors[CCs_new[i].nodes()] = 0

facets_gray = np.where(colors == 0)[0]
facets_red  = np.where(colors == 1)[0]

F_colored = np.hstack((F, np.ones((F.shape[0], 4))))
F_colored[facets_gray, 4:7] = 200*F_colored[facets_gray, 4:7]
F_colored[facets_red, 4]    = 255*F_colored[facets_red, 4]
F_colored[facets_red, 5]    = 0
F_colored[facets_red, 6]    = 0
F_colored[:, 7] = 255

# Final projection:
with open(output_folder + "projected_channel.off", 'w') as output_mesh:
    output_mesh.write("COFF\n\n")
    output_mesh.write(f"{V.shape[0]} {F.shape[0]} 0\n")
    
    for v in V:
        output_mesh.write(f"{v[0]} {v[1]} {v[2]}\n")
        
    for f in F_colored:
        output_mesh.write(f"{int(f[0])} {int(f[1])} {int(f[2])} {int(f[3])} {int(f[4])} {int(f[5])} {int(f[6])} {int(f[7])}\n")



# CHANNEL EXTRACTION:
F_red = F_colored[facets_red, 1:4]
V_red = V[np.unique(F_red.flatten()).astype(int).tolist(), :]

with open(output_folder + "channel.off", 'w') as output_mesh:

    knn = NearestNeighbors(n_neighbors=1).fit(V_red)
    _, idx_1 = knn.kneighbors(V[F_red[:, 0].astype(int).tolist(), :])
    _, idx_2 = knn.kneighbors(V[F_red[:, 1].astype(int).tolist(), :])
    _, idx_3 = knn.kneighbors(V[F_red[:, 2].astype(int).tolist(), :])
    idx_1 = np.squeeze(np.array(idx_1)).tolist()
    idx_2 = np.squeeze(np.array(idx_2)).tolist()
    idx_3 = np.squeeze(np.array(idx_3)).tolist()

    output_mesh.write("OFF\n\n")
    output_mesh.write(f"{V_red.shape[0]} {len(idx_1)} 0\n")
    
    for v in V_red:
        output_mesh.write(f"{v[0]} {v[1]} {v[2]}\n")

    for i1, i2, i3 in zip(idx_1, idx_2, idx_3):
        output_mesh.write(f"3 {i2} {i1} {i3}\n")


end = time.time()
print(f"Elapsed time (in seconds): {end - start}")
