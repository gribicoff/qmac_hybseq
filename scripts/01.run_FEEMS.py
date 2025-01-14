## IMPORT LIBRARIES
# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import os

# viz
import matplotlib.pyplot as plt
from matplotlib import gridspec
import cartopy.crs as ccrs
import geopandas as gpd

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"

# Change your directory
os.chdir("/Users/kieranalthaus/Dropbox/Mac/Desktop/Projects/QMAC_FEEMS")

## READ DATA
# Read in PLINK formatted data. NOTE: FEEMS doesn't like missing SNPs,
# so here you should compute the mean SNP value at all missing locations.
# Further, FEEMS doesn't like monomorphic sites, so get rid of those too

# Path to data
path_to_plink = "DATA/new_data/"
# Input the file prefix and load files
(bim, fam, G) = read_plink("{}/qmac_feems".format(path_to_plink)) # 'qmac_feems' is the title of my files. Input your own

# mean impute missing data
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

# setup graph
coord = np.loadtxt("{}/qmac.coord".format(path_to_plink))  # sample coordinates. 'qmac.coord' is my file prefix
outer = np.loadtxt("{}/qmac.outer".format(path_to_plink))  # outer coordinates. 'qmac.outer' is my file prefix
grid_path = "{}/grid_100.shp".format(path_to_plink)  # path to discrete global grid

# graph input files
outer, edges, grid, _ = prepare_graph_inputs(coord=coord, 
                                             ggrid=grid_path,
                                             translated=True, 
                                             buffer=0,
                                             outer=outer)

# Set up sptaial graph
sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)
projection = ccrs.EquidistantConic(central_longitude=-108.842926, central_latitude=66.037547)

## Conduct Cross-Validation to find optimal lambda
# Define range of lambda values
lambda_grid = np.geomspace(10e-5, 10, 20)[::-1]

# Run cross validation
cv_err = run_cv(sp_graph, lambda_grid,n_folds=sp_graph.n_observed_nodes,factr=1e10)

# Average over folds
mean_cv_err = np.mean(cv_err, axis=0)

# ARGMIN of cv error
lambd_cv = float(lambda_grid[np.argmin(mean_cv_err)])

# Plot and save cross-validation error
fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lambda_grid), mean_cv_err, ".");
ax.set_xlabel("log10(lambda)");
ax.set_ylabel("L2 CV Error");
fig.text(np.log10(lambd_cv), 0, "0.736", transform=ax.transAxes);
ax.axvline(np.log10(lambd_cv), color = "orange")
plt.savefig("OUT/20250113_qmac_cv_error.png")

## Run FEEMS 
fig = plt.figure(constrained_layout=True, dpi=300, figsize=(6, 6))
spec = plt.GridSpec(ncols=1, nrows=2, figure=fig)
lamb = float(lambd_cv)
sp_graph.fit(lamb = lamb)

## EXPORT DATA
# Set up edge and edge weight dfs
edges_df = pd.DataFrame(sp_graph.edges, columns = ['node1', 'node2'])
weights_df = pd.DataFrame(sp_graph.w, columns = ['weight'])
# Save edge weights
edges_weights_df = pd.concat([edges_df, weights_df], axis = 1)

# Set up sample positions
sample_assignments = pd.DataFrame(
    data = sp_graph.n_samples_per_obs_node_permuted,
    columns = ['n_samples']
    )
sample_pos = pd.DataFrame(sp_graph.sample_pos, columns = ['longitude', 'latitude'])
samples_df = pd.concat([sample_pos, sample_assignments], axis = 1)
# Save sample positions
samples_df.to_csv('OUT/sample_positions.csv', index = False)

# Set up node positions
node_pos_df = pd.DataFrame(sp_graph.node_pos, columns = ['longitude', 'latitude'])
node_pos_df.to_csv('OUT/node_positions.csv', index = False)