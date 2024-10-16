import numpy as np
import pandas as pd
import matplotlib.pyplot
import matplotlib.pyplot as plt
import sys
import os
import dynamo as dyn
import seaborn as sns
import anndata as ann
adata = dyn.read_h5ad("/data1/fanlab/wry/final/n00h/n00htotal.h5ad")
geneset = pd.read_csv("/data1/fanlab/wry/final/n00h/qcfil.csv",index_col=1)
geneset = list(geneset.index)
geneset = list(geneset)
print(geneset)
dyn.pp.recipe_monocle(adata,tkey="label_time",normalized=True,experiment_type="one-shot",genes_to_use=geneset,n_top_genes=len(geneset),maintain_n_top_genes=True)
#dyn.tl.reduceDimension(adata)
dyn.tl.dynamics(adata)
n=pd.DataFrame(adata.var)
n.to_csv("/data1/fanlab/wry/final/n00h/var_n00h_1.csv",index=True,sep=",")
