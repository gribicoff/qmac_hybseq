#!/usr/bin/env python
# coding: utf-8

# ### Import libraries
# 
# Make sure you're using the correct conda environment

# In[37]:


# base
import numpy as np
import pkg_resources
from sklearn.impute import SimpleImputer
from pandas_plink import read_plink
import os
import statsmodels.api as sm

# viz
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib import gridspec
import geopandas as gpd

# feems
from feems.utils import prepare_graph_inputs
from feems import SpatialGraph, Viz
from feems.cross_validation import run_cv

# change matplotlib fonts
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.sans-serif"] = "Arial"


# In[ ]:

# #### Read Data

# Read the PLINK formatted data
# NOTE: FEEMS don't like missing SNPSs, so here you'll compute the mean SNP value at all missing locations
# ALSO, feems doesn't want monomorphic sites. If you haven't already done so, use the "--maf" command in plink to remove those sites

# In[29]:


# Path to data
path_to_plink = "./"
# Input the file prefix and load files
(bim, fam, G) = read_plink("{}/qmac_feems".format(path_to_plink)) # 'qmac_maff_feems' is the title of my files. Input your own

# mean impute missing data
imp = SimpleImputer(missing_values=np.nan, strategy="mean")
genotypes = imp.fit_transform((np.array(G)).T)

# Print
print("n_samples={}, n_snps={}".format(genotypes.shape[0], genotypes.shape[1]))


# For the following step, you'll need to have downloaded some grid shapefile as provided by the github page "https://github.com/NovembreLab/feems/tree/main/feems/data". I use the grid_100.shp below.

# In[30]:


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


# In[31]:


sp_graph = SpatialGraph(genotypes, coord, grid, edges, scale_snps=True)


# In[90]:


projection = ccrs.EquidistantConic(central_longitude=-95.0, central_latitude=40.0)


# Draw out the map with grids and points plotted.

# In[50]:


fig = plt.figure(dpi = 100)
ax = fig.add_subplot(1,1,1, projection=projection)
v = Viz(ax, sp_graph, projection=projection, edge_width=.5,
        edge_alpha=1, sample_pt_size=10, edge_zorder=100,
        obs_node_size=7.5, sample_pt_color="black",
        cbar_font_size=10)
v.draw_map()
v.ax.add_feature(cfeature.STATES, linewidth=0.25, edgecolor="#636363", zorder=0)
v.draw_samples()
v.draw_edges(use_weights=False)
v.draw_obs_nodes(use_ids=False)


# Use LOO cross-validation to determine the best value of lambda (smoothing parameter):

# In[51]:


lamb_grid = np.geomspace(1e-6, 1e2, 20)[::-1]
cv_err = run_cv(sp_graph, lamb_grid, n_folds=sp_graph.n_observed_nodes, factr=1e10)
mean_cv_err = np.mean(cv_err, axis=0)
lamb_cv = float(lamb_grid[np.argmin(mean_cv_err)])


# Plot CV error as a function of lambda:

# In[86]:


fig, ax = plt.subplots(dpi=300)
ax.plot(np.log10(lamb_grid), mean_cv_err, ".");
ax.set_xlabel("$\mathrm{log}_{10}(\lambda)$");
ax.set_ylabel("LOO CV Error");
ax.axvline(np.log10(lamb_cv), color = "orange", linestyle="--")
ax.text(0.03, 0.03, "$\lambda_{{err_{{min}}}}=${}".format(lamb_cv), transform=ax.transAxes)


# Save figure:

# In[87]:


fig.savefig(fname= "qmac_feems_cverror.pdf")


# Generate FEEMS output with best lambda value, as chosen by CV:

# In[88]:


fig = plt.figure(constrained_layout=True, dpi=300, figsize=(6, 6))
spec = plt.GridSpec(ncols=1, nrows=2, figure=fig)
sp_graph.fit(lamb = lamb_cv)


# Now plot the result

# In[109]:


fig = plt.figure(dpi=300)
ax = fig.add_subplot(1,1,1,projection=projection)
v = Viz(ax, sp_graph, projection = projection, edge_width=.5,
        edge_alpha = 1, edge_zorder=100, sample_pt_size = 20,
        obs_node_size = 7.5, sample_pt_color = "black",
        cbar_font_size = 7, cbar_ticklabelsize = 7, cbar_loc = "lower right")
v.draw_map()
v.ax.add_feature(cfeature.STATES, linewidth=0.25, edgecolor="#636363", zorder=0)
v.draw_edges(use_weights=True)
v.draw_obs_nodes(use_ids=False)
v.draw_edge_colorbar()



# Save as PDF:

# In[110]:


fig.savefig(fname= "qmac_feems.pdf")

