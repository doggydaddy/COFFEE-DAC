# Investigating edge bundling

# imports


# general
import numpy as np

# data loading
import pandas as pd
import nibabel as nib

# clustering libraries
import scipy.cluster.hierarchy
import sklearn.cluster
from sklearn.cluster import DBSCAN
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

# distance eval
import matplotlib.pyplot as plt

# visualization
import pyvista as pv
import matplotlib as mpl

# misc
from tqdm import tqdm

# -----------------
# support functions

def euc_dist(c1, c2):
    '''
    Given 2 points, find eucledian distance between them.
    '''
    
    dx = pow( abs(c1[0] - c2[0]) , 2)
    dy = pow( abs(c1[1] - c2[1]) , 2)
    dz = pow( abs(c1[2] - c2[2]) , 2)
    
    return( np.sqrt(dx+dy+dz) )  

def edge_dist(edge_a, edge_b, flag, directional=False):

    edge_a_ep1 = edge_a[0:3]
    edge_a_ep2 = edge_a[3:6]
    edge_b_ep1 = edge_b[0:3]
    edge_b_ep2 = edge_b[3:6]

    if directional:
        d1 = euc_dist(edge_a_ep1, edge_b_ep1)
        d2 = euc_dist(edge_a_ep1, edge_b_ep2)
        d3 = euc_dist(edge_a_ep2, edge_b_ep1)
        d4 = euc_dist(edge_a_ep2, edge_b_ep2)
        dists = (d1, d2, d3, d4)
    else:
        d1 = euc_dist(edge_a_ep1, edge_b_ep1)
        d4 = euc_dist(edge_a_ep2, edge_b_ep2)
        dists = (d1, d4)
         
    if flag == 'min':
        outp = np.min(dists)
    elif flag == 'max':
        outp = np.max(dists)
    elif flag == 'mean':
        outp = np.mean(dists)
    else:
        print("edge_dist error: flag must be either 'min', 'max', or 'mean'")
        return( -1 )
    
    return outp



# --------------------------------------------
# edge-bundling: First hierarchical clustering


# calculate distance matrix for first hclust
def h1_dist(edges, flag):

    # total number of edges
    N = edges.shape[0]
    
    # allocate distance matrix
    dist = np.zeros((N, N))

    # foreach pair of edges (a, b):
    for a in range(0, N, 1):
        for b in range(a, N, 1):
            
            edge_a = edges[a, :]
            edge_b = edges[b, :]
            
            dist[a, b] = edge_dist(edge_a, edge_b, flag, False)
            
            if (a != b ):
                dist[b, a] = dist[a, b]

    return( dist )


# perform first hclust
def hc1(edges, dist, flag, nr_clusters):
    ncl = nr_clusters 
    if flag == 'complete':
        clustering = AgglomerativeClustering(n_clusters=ncl, affinity='precomputed', linkage='complete').fit(dist)
    elif flag == 'average':
        clustering = AgglomerativeClustering(n_clusters=ncl, affinity='precomputed', linkage='average').fit(dist)
    else:
        print("Unknown flag, must be either 'complete' or 'average'")
        return( -1 )
    
    # loading clustering into edge list
    edges_out = edges
    edges_out[:, 6] = clustering.labels_
    return( edges_out )


# analytics for first hclust: calculate max-min metric
def maxmin(edges, cluster_nr):
    
    edges_in_cl = edges[edges[:,6]==cluster_nr, :]
    min_dists = h1_dist(edges_in_cl, 'min')
    max_min = np.max(min_dists)
    return(max_min)



# ---------------------------------------------
# functional connectivity network construction:
# second hierarchical clustering


def h2_dist(edges, idx1, idx2):
    
    edges_1 = edges[edges[:,6]==idx1, :]
    edges_2 = edges[edges[:,6]==idx2, :]
    
    ne1 = edges_1.shape[0]
    ne2 = edges_2.shape[0]
    
    dmins = np.zeros(ne1*ne2)
    for i in range(ne1):
        for j in range(ne2):
            
            edge_i = edges_1[i, :]
            edge_j = edges_2[j, :]
            
            dmins[i*ne2+j] = edge_dist(edge_i, edge_j, 'min', False)
            
    return( np.min(dmins) )


def hc2(edges, dist, flag, nr_networks):
    
    if flag == 'single':
        clustering = AgglomerativeClustering(n_clusters=nr_networks, affinity='precomputed', linkage='single').fit(dist)
    elif flag == 'average':
        clustering = AgglomerativeClustering(n_clusters=nr_networks, affinity='precomputed', linkage='average').fit(dist)
    elif flag == 'complete':
        clustering = AgglomerativeClustering(n_clusters=nr_networks, affinity='precomputed', linkage='complete').fit(dist)
    else:
        print("Unknown flag, must be either 'single' or 'average'")
        return( -1 )
    
    # put network information to edge list
    N = edges.shape[0]
    network_labels = np.zeros(N)
    for i in range(N):
        network_labels[i] = clustering.labels_[ edges[i, 6] ]

    edges_net = np.c_[edges, network_labels]
    return(edges_net)



# -------------
# Visualization


def plot_edge(plotter, edge):
    this_color = my_colormap(int(edge[6]))
    this_opa = 1
    
    # endpoints for the spline
    a = edge[0:3]
    b = edge[3:6]
    
    # constructing midpoint with offset as "midpoint" to make the splines a bit nicer.
    offset_i = np.abs(a[0] - b[0]) * 0.1
    offset_j = np.abs(a[1] - b[1]) * 0.1
    offset_k = np.abs(a[2] - b[2]) * 0.1
    m = [((a[0]+b[0])/2)+offset_i, 
         ((a[1]+b[1])/2)+offset_j, 
         ((a[2]+b[2])/2)+offset_k]
    
    # defining the spline
    spline_points = np.row_stack((a, m, b))
    
    N = 32 # how many spline segments
    spline = pv.Spline( spline_points, N )
    spline['scalars'] = np.full(N, int(edge[6]))
    
    # adding endpoints
    d=np.sqrt(3)/4
    a_box = pv.Box([a[0]-d, a[0]+d, a[1]-d, a[1]+d, a[2]-d, a[2]+d])
    b_box = pv.Box([b[0]-d, b[0]+d, b[1]-d, b[1]+d, b[2]-d, b[2]+d])
    plotter.add_mesh(a_box, color=this_color)
    plotter.add_mesh(b_box, color=this_color)

    plotter.add_mesh(spline, color=this_color, render_lines_as_tubes=True, line_width=3, opacity=1.0, show_scalar_bar=True)
    

def plot_net(plotter, edge):
    this_color = my_colormap(int(edge[7]))
    this_opa = 1
    
    # endpoints for the spline
    a = edge[0:3]
    b = edge[3:6]
    
    # constructing midpoint with offset as "midpoint" to make the splines a bit nicer.
    offset_i = np.abs(a[0] - b[0]) * 0.1
    offset_j = np.abs(a[1] - b[1]) * 0.1
    offset_k = np.abs(a[2] - b[2]) * 0.1
    m = [((a[0]+b[0])/2)+offset_i, 
         ((a[1]+b[1])/2)+offset_j, 
         ((a[2]+b[2])/2)+offset_k]
    
    # defining the spline
    spline_points = np.row_stack((a, m, b))
    
    # adding endpoints
    d=np.sqrt(3)/4
    a_box = pv.Box([a[0]-d, a[0]+d, a[1]-d, a[1]+d, a[2]-d, a[2]+d])
    b_box = pv.Box([b[0]-d, b[0]+d, b[1]-d, b[1]+d, b[2]-d, b[2]+d])
    plotter.add_mesh(a_box, color=this_color)
    plotter.add_mesh(b_box, color=this_color)
    
    N = 32 # how many spline segments
    spline = pv.Spline( spline_points, N )
    spline['scalars'] = np.full(N, int(edge[7]))

    plotter.add_mesh(spline, color=this_color, render_lines_as_tubes=True, line_width=3, opacity=0.3, show_scalar_bar=True)
    
    
    
# Define the colors we want to use, currently supports up to 13 * 2 = 26 colours
# For more colors, https://www.rapidtables.com/web/color/RGB_Color.html

from matplotlib.colors import ListedColormap

one = np.array([153/256, 0.0, 0.0, 1.0])
two = np.array([153/256, 76/256, 0.0, 1.0])
three = np.array([153/256, 153/256, 0.0, 1.0])
four = np.array([76/256, 153/256, 0.0, 1.0])
five = np.array([0/256, 153/256, 0.0, 1.0])
six = np.array([0/256, 153/256, 76/256, 1.0])
seven = np.array([0/256, 153/256, 153/256, 1.0])
eight = np.array([0/256, 76/256, 153/256, 1.0])
nine = np.array([0/256, 0/256, 153/256, 1.0])
ten = np.array([76/256, 0/256, 153/256, 1.0])
eleven = np.array([153/256, 0/256, 153/256, 1.0])
twelve = np.array([153/256, 0/256, 76/256, 1.0])
thirteen = np.array([64/256, 64/256, 164/256, 1.0])

fourteen = np.array([255/256, 102/256, 102/256, 1.0])
fifteen = np.array([255/256, 178/256, 102/256, 1.0])
sixteen = np.array([255/256, 255/256, 102/256, 1.0])
seventeen = np.array([178/256, 255/256, 102/256, 1.0])
eighteen = np.array([102/256, 255/256, 102/256, 1.0])
nineteen = np.array([102/256, 255/256, 178/256, 1.0])
twenty = np.array([102/256, 255/256, 255/256, 1.0])
twentyone = np.array([102/256, 178/256, 255/256, 1.0])
twentytwo = np.array([102/256, 102/256, 255/256, 1.0])
twentythree = np.array([178/256, 102/256, 255/256, 1.0])
twentyfour = np.array([255/256, 102/256, 255/256, 1.0])
twentyfive = np.array([255/256, 102/256, 178/256, 1.0])
twentysix = np.array([192/256, 192/256, 192/256, 1.0])

twentyseven = np.array([255/256, 0/256, 0/256, 1.0])
twentyeight = np.array([255/256, 128/256, 0/256, 1.0])
twentynine = np.array([255/256, 255/256, 0/256, 1.0])
thirty = np.array([128/256, 255/256, 0/256, 1.0])
thirtyone = np.array([0/256, 255/256, 0/256, 1.0])
thirtytwo = np.array([0/256, 255/256, 128/256, 1.0])
thirtythree = np.array([0/256, 255/256, 255/256, 1.0])
thirtyfour = np.array([0/256, 128/256, 255/256, 1.0])
thirtyfive = np.array([0/256, 0/256, 255/256, 1.0])
thirtysix = np.array([128/256, 0/256, 255/256, 1.0])
thirtyseven = np.array([255/256, 0/256, 255/256, 1.0])
thirtyeight = np.array([255/256, 0/256, 128/256, 1.0])
thirtynine = np.array([128/256, 128/256, 128/256, 1.0])



mapping = np.linspace(1, 14, 39)
newcolors = np.empty((39, 4))
newcolors[0] = one
newcolors[1] = two
newcolors[2] = three
newcolors[3] = four
newcolors[4] = five
newcolors[5] = six
newcolors[6] = seven
newcolors[7] = eight
newcolors[8] = nine
newcolors[9] = ten
newcolors[10] = eleven
newcolors[11] = twelve
newcolors[12] = thirteen
newcolors[13] = fourteen
newcolors[14] = fifteen
newcolors[15] = sixteen
newcolors[16] = seventeen
newcolors[17] = eighteen
newcolors[18] = nineteen
newcolors[19] = twenty
newcolors[20] = twentyone
newcolors[21] = twentytwo
newcolors[22] = twentythree
newcolors[23] = twentyfour
newcolors[24] = twentyfive
newcolors[25] = twentysix
newcolors[26] = twentyseven
newcolors[27] = twentyeight
newcolors[28] = twentynine
newcolors[29] = thirty
newcolors[30] = thirtyone
newcolors[31] = thirtytwo
newcolors[32] = thirtythree
newcolors[33] = thirtyfour
newcolors[34] = thirtyfive
newcolors[35] = thirtysix
newcolors[36] = thirtyseven
newcolors[37] = thirtyeight
newcolors[38] = thirtynine

# Make the colormap from the listed colors
my_colormap = ListedColormap(newcolors)



# Investigating edge bundling
NR_BUNDLES = 44
NR_NETWORKS = 5

edges_csv_file = 'edges_to_plot_ijk.csv'
edges_ijk = pd.read_csv(edges_csv_file)
edges = edges_ijk.to_numpy()

ep_max = h1_dist(edges, "max") # default
edges = hc1(edges, ep_max, 'complete', NR_BUNDLES)

# calculate distance between edge bundles
# Closest distance between any edge in bundle A any edge in bundle B
dist2 = np.zeros([NR_BUNDLES, NR_BUNDLES])
for i in range(NR_BUNDLES):
    for j in range(i, NR_BUNDLES):
        if i == j:
            dist2[i, j] = 0.
        else:
            dist2[i, j] = h2_dist(edges, i, j)
            dist2[j, i] = dist2[i, j]

edges_net = hc2(edges, dist2, 'average', NR_NETWORKS)



# visualize stl
p = pv.Plotter()
reader = pv.get_reader("test.stl")
mesh = reader.read()
p.add_mesh(mesh, show_edges=True, edge_color="gray", color='white', point_size=3, render_points_as_spheres=True, opacity=0.1, show_scalar_bar=True)


# plot an network 
'''
network_id = 4
edges_to_plot = edges_net[edges_net[:, 7] == network_id, :]
ne = edges_to_plot.shape[0]
for i in range(ne):
    plot_net(p, edges_to_plot[i,:])
'''
    
# plot networks 
network_ids = [2, 3]
#network_ids = range(3, 20)
for i in network_ids:
    edges_to_plot = edges_net[edges_net[:, 7] == i, :]
    ne = edges_to_plot.shape[0]
    for j in range(ne):
        plot_net(p, edges_to_plot[j,:])
    
'''
# plot all networks
edges_to_plot = edges_net[edges_net[:, 7] >= 0, :]
ne = edges_to_plot.shape[0]
for i in range(ne):
    plot_net(p, edges_to_plot[i,:])
'''   
  
p.show(jupyter_backend='trame')