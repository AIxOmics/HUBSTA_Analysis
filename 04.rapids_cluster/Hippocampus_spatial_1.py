import numpy as np
import scanpy as sc
import rapids_singlecell as rsc
import scanpy as sc
import rapids_singlecell as rsc
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar 
from matplotlib.colors import ListedColormap, rgb2hex
import pickle
import numpy as np
import warnings
import pandas as pd
warnings.filterwarnings('ignore')


adata = sc.read_h5ad('/data/work/05.cluster/FuseMap/0106/Hippocampus.h5ad')


adata.obsm['X_dmt_highdim'] = np.load("/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/X_dmt_highdim.npy")
adata.obsm['X_dmt'] = np.load("/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/X_dmt.npy")
rsc.pp.neighbors(adata,use_rep='X_dmt_highdim',key_added='dmt_nn') # 50维
rsc.tl.leiden(adata, neighbors_key='dmt_nn',key_added='dmt_leiden')


number_of_colors = len(set(adata.obs['dmt_leiden']))
colors = np.random.rand(number_of_colors, 3)
custom_cmap = ListedColormap(colors)
hex_colors = [rgb2hex(color) for color in custom_cmap.colors]
colormap = {str(i):hex_colors[i] for i in range(number_of_colors) }

with open('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/colormap20250108_1.pkl', 'wb') as f:
    pickle.dump(colormap, f)

with open('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/colormap20250108_1.pkl', 'rb') as f:
    colormap = pickle.load(f)

sc.pl.embedding(adata, basis = 'X_dmt', color = 'slice_code', show = False)
plt.savefig('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/batch_20250108_1.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()


sc.pl.embedding(adata, basis = 'X_dmt', color = 'dmt_leiden',palette=colormap, show = False)
plt.savefig('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/dmt_leiden_20250108_1.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()

names = [
    '12_B03605F3G5_WT202403310048.h5ad',
    '13_B03612A1C3_WT202403310056.h5ad',
    '14_A03591A1C3_WT202403310045.h5ad',
    '16_A03592A4C6_WT202403310044.h5ad',
    '18_B03602C4D6_WT202405020031.h5ad',
    '20_B03606F3G5_WT202405020032.h5ad',
    '22_B03606C4E6_WT202403310050.h5ad',
    '23_B03609A4D6_WT202404150263.h5ad',
    '27_B03610C1E3_WT202403310051.h5ad',
    '31_B03619A1D3_WT202403310052.h5ad',
    '35_B03619E4G6_WT202403310053.h5ad',
    '39_A03589A1D4_WT202403310046.h5ad',
    '43_A03590E1G4_WT202403310064.h5ad',
    '47_A03593C1F3_WT202403310068.h5ad',
    'B03607C4E6_WT2024071214941.h5ad',
]



fig = plt.figure(figsize=(32, 12))  

gs = GridSpec(2, 8, figure=fig)  
count = 0
for name in names:
    adata_temp = adata[adata.obs['slice_code'] == name].copy()
    adata_temp.obsm['align_spatial_2d'][:, 1] = -adata_temp.obsm['align_spatial_2d'][:, 1]

    row = (count // 8) + 1  
    col = count % 8      
    count +=1
    ax = fig.add_subplot(gs[row-1, col])
    sc.pl.embedding(adata_temp, basis="align_spatial_2d", color='dmt_leiden', 
                    show=False, s=1, palette=colormap, title='', legend_loc=None, ax=ax)
    ax.axis('off')
    ax.set_aspect('equal')

    scalebar = ScaleBar(0.0097, "mm", fixed_value=1, location = 'lower left', frameon = False,);
    ax.add_artist(scalebar)

plt.tight_layout()  

plt.savefig('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/dmt_leiden_all_20250108_1.png', dpi = 450, bbox_inches = 'tight')
plt.close()

adata.write('/data/work/05.cluster/FuseMap/0106/Hippocampus_latent_embeddings_all_spatial_pretrain/dmt_leiden_20250108_1.h5ad')