#!/usr/bin/env python
# coding: utf-8

# In[1]:


import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import numba

print(scv.__version__)
print(np.__version__)
print(pd.__version__)
print(sc.__version__)
print(numba.__version__)


# In[2]:


adata1 = scv.read("/banach1/rq25/elsie_velocity/Old_Females.h5ad")
adata2 = scv.read("/banach1/rq25/elsie_velocity/Old_Males.h5ad")
adata3 = scv.read("/banach1/rq25/elsie_velocity/Young_Females.h5ad")
adata4 = scv.read("/banach1/rq25/elsie_velocity/Young_Males.h5ad")

adata = adata1.concatenate(adata2)
adata = adata.concatenate(adata3)
adata = adata.concatenate(adata4)
adata


# In[3]:


adata.obs["seurat_clusters"].value_counts()


# In[4]:


adata = adata[adata.obs["seurat_clusters"]!= 9]
adata


# In[5]:


scv.pl.proportions(adata)


# In[6]:


adata.obs["seurat_clusters"].value_counts()


# In[7]:


scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=5000)
scv.pp.log1p(adata)


# In[8]:


scv.pp.moments(adata, n_pcs=30, n_neighbors=30)


# ## Steady-state Model
# - https://www.nature.com/articles/s41586-018-0414-6

# In[9]:


scv.tl.velocity(adata)


# In[10]:


scv.tl.velocity_graph(adata)


# In[11]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters", dpi = 200, legend_loc = "right margin", save = "../Plot/5.RNA_velocity/steady_state_seurat_cluster.png", show = False)
scv.pl.velocity_embedding_stream(adata, basis='umap', color="orig.ident", dpi = 200, legend_loc = "right margin",save = "../Plot/5.RNA_velocity/steady_state_orig_ident.png", show = False)


# In[12]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters", dpi = 200, legend_loc = "right margin")


# In[13]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="orig.ident", dpi = 200, legend_loc = "right margin")


# In[14]:


scv.tl.rank_velocity_genes(adata, groupby='seurat_clusters', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head(20)


# In[15]:


scv.pl.velocity(adata, ['Itga4'], ncols=1, color = 'seurat_clusters', basis = "umap", dpi = 100)


# In[16]:


for i in [str(x) for x in np.unique(adata.obs['seurat_clusters'])]:
    scv.pl.velocity(adata, df[i][:20], ncols=2, color = 'seurat_clusters', basis = "umap", dpi = 100, save = "../Plot/5.RNA_velocity/steady_state_model/steady_state_model_cluster" + i + "_top20.png", show = False)


# In[17]:


scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], basis = "umap")


# In[18]:


scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', basis = "umap")


# In[22]:


scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], basis = "umap",save = "../Plot/5.RNA_velocity/steady_state_velocity_length_confidence.png")
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', basis = "umap",save = "../Plot/5.RNA_velocity/steady_state_velocity_pseudotime.png")


# ## Dynamic Model
# - (https://www.nature.com/articles/s41587-020-0591-3)

# In[23]:


scv.tl.recover_dynamics(adata)


# In[24]:


scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, basis = "umap")


# In[25]:


scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)


# In[26]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters", dpi = 200, legend_loc = "right margin", save = "../Plot/5.RNA_velocity/dynamic_seurat_cluster.png", show = False)
scv.pl.velocity_embedding_stream(adata, basis='umap', color="orig.ident", dpi = 200, legend_loc = "right margin",save = "../Plot/5.RNA_velocity/dynamic_orig_ident.png", show = False)


# In[27]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="seurat_clusters", dpi = 200, legend_loc = "right margin")


# In[28]:


scv.pl.velocity_embedding_stream(adata, basis='umap', color="orig.ident", dpi = 200, legend_loc = "right margin")


# In[29]:


scv.tl.rank_dynamical_genes(adata, groupby='seurat_clusters')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(20)


# In[30]:


scv.pl.velocity(adata, df['4'][:10], ncols=2, color = 'seurat_clusters', basis = "umap", dpi = 100)


# In[31]:


for i in [str(x) for x in np.unique(adata.obs['seurat_clusters'])]:
    scv.pl.velocity(adata, df[i][:20], ncols=2, color = 'seurat_clusters', basis = "umap", dpi = 100, save = "../Plot/5.RNA_velocity/dynamic_model/dynamic_model_cluster" + i + "_top20.png", show = False)


# In[33]:


scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, basis = "umap",save = "../Plot/5.RNA_velocity/latent_time.png")


# In[ ]:




