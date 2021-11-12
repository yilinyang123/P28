import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

@st.cache(allow_output_mutation=True, show_spinner=False)
def load_data():
	adata = sc.read("adata_P28.h5ad")
	#sc.pp.normalize_total(adata) #Normalize data by the median counts per single cell library
	#adata.X = np.arcsinh(adata.X).copy() #Transform the data by using an inverse hyperbolic sine transform, this eliminates the need for adding a pseudocount 
	#sc.pp.scale(adata) #Scale and center the data for interpretability
	return adata
adata = load_data()

gene_list = list(adata.var_names)
st.set_option('deprecation.showPyplotGlobalUse', False)
st.title('thymus')

selectbox = st.sidebar.selectbox(
    'Section',
    ('UMAP', 'Genotype')
)

user_input = st.sidebar.multiselect("gene",gene_list,['Cd83'])
if len(user_input) >0:
	user_input = user_input[0]
else: user_input = 'Cd83'
coord = adata.obsm['X_umap']
a = adata[: , user_input].X.toarray().astype('float')
a_min = np.round_(np.min(a), decimals=1)
a_max = np.round_(np.max(a), decimals=1)
x = st.slider("Scale", 0.0, 3.0, (float(a_min), float(a_max)), 0.1)
plt.rcParams.update({
    "figure.facecolor":  (0.0, 0.0, 0.0, 0.0),  # red   with alpha = 30%
    #"axes.facecolor":    (0.0, 0.0, 0.0, 0.0),  # green with alpha = 50%
})

@st.cache(allow_output_mutation=True, show_spinner=False)
def get_fig():
	fig,ax = plt.subplots(figsize=(15,6))
	sc = ax.scatter(coord[:, 0], coord[:, 1], c = a, cmap = 'magma', s = 3,alpha=0.65, vmin = x[0], vmax = x[1])
	cb = fig.colorbar(sc, ax=ax)
	ax.axis('off')
	ax.set_title(user_input,color='white')
	cb.ax.yaxis.set_tick_params(color='white')
	plt.setp(plt.getp(cb.ax.axes, 'yticklabels'), color='white')
	return fig
	# fig = sc.pl.umap(adata, color=["leiden",user_input],
 #           use_raw=False,wspace=.1, hspace=0.2,
 #           frameon=False,add_outline=True,
 #           legend_loc='on data',
 #           title=["Clusters",user_input],
 #           palette='tab20')
	# return fig

@st.cache(allow_output_mutation=True, show_spinner=False)
def get_fig2():
	color_map = {0: 'r', 1:'b'}
	colors = adata.obs['genotype_new'].map(color_map)
	fig,ax = plt.subplots(figsize=(15,6))
	ax.scatter(coord[:, 0], coord[:, 1], c = colors, cmap = 'viridis', s = 2,alpha=0.65,label=adata.obs['genotype_new'])
	#fig.colorbar(sc, ax=ax2)
	ax.axis('off')
	ax.set_title('Genotype')
	return fig

if selectbox == 'UMAP':
	st.pyplot(get_fig())

if selectbox == 'Genotype':
	st.pyplot(get_fig2())