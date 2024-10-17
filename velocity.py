import numpy as np
import pandas as pd
import matplotlib.pyplot
import matplotlib.pyplot as plt
import sys
import os
import dynamo as dyn
import seaborn as sns
import anndata as ann
adata = dyn.read_h5ad("/data1/fanlab/wry/final/B1/0h/Btotal0h.h5ad")
dyn.pp.recipe_monocle(adata,tkey="label_time",normalized=True,experiment_type="one-shot",n_top_genes=1000)
dyn.tl.reduceDimension(adata)
dyn.tl.dynamics(adata, one_shot_method="sci_fate", model="stochastic")
dyn.tl.cell_velocities(adata,basis="umap_ori")
color_list = {"B":"#81b8df",
              "Plasma":"#fe817d"}
dyn.pl.streamline_plot(adata,figsize=(7, 4),color_key=color_list,color="cell_type",inverse=True,show_legend=False,basis="umap_ori",save_show_or_return="return")
dyn.pl.save_fig(path="/data1/fanlab/wry/final/B1/0h/",prefix="B_0h_umap",ext="svg",close=True,verbose=True)
