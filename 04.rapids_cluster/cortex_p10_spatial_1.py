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


adata = sc.read_h5ad('/data/work/05.cluster/FuseMap/1230/cortex_p10.h5ad')


adata.obsm['X_dmt_highdim'] = np.load("/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/X_dmt_highdim.npy")
adata.obsm['X_dmt'] = np.load("/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/X_dmt.npy")
adata = adata[adata.obs['slice_code'] != 'A03591D4E5_WT2024071215074.h5ad']
rsc.pp.neighbors(adata,use_rep='X_dmt_highdim',key_added='dmt_nn') # 50维
rsc.tl.leiden(adata, neighbors_key='dmt_nn',key_added='dmt_leiden')


number_of_colors = len(set(adata.obs['dmt_leiden']))
colors = np.random.rand(number_of_colors, 3)
custom_cmap = ListedColormap(colors)
hex_colors = [rgb2hex(color) for color in custom_cmap.colors]
colormap = {str(i):hex_colors[i] for i in range(number_of_colors) }

with open('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/colormap20250108_1.pkl', 'wb') as f:
    pickle.dump(colormap, f)

with open('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/colormap20250108_1.pkl', 'rb') as f:
    colormap = pickle.load(f)

sc.pl.embedding(adata, basis = 'X_dmt', color = 'slice_code', show = False)
plt.savefig('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/batch_20250108_1.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()


sc.pl.embedding(adata, basis = 'X_dmt', color = 'dmt_leiden',palette=colormap, show = False)
plt.savefig('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/dmt_leiden_20250108_1.png', dpi=300, bbox_inches='tight', pad_inches=0)
plt.close()

names = [
 # '1_C03627A3_WT202403110495.h5ad',
 '2_C03627A6_WT202403110494.h5ad',
 '3_C03627F1_WT202403110496.h5ad',
 '4_C03627G1_WT202403110489.h5ad',
 '5_D03657E5_WT202403110487.h5ad',
 '6_C02941A1B2_WT202403310083.h5ad',
 '7_D03473F1G2_WT202403310080.h5ad',
 '8_B03422F5G6_WT202403310079.h5ad',
 '9_D03466F1G2_WT202403310042.h5ad',
 '10_D03470F1G2_WT202403310077.h5ad',
 '11_B03612A4C6_WT202403310047.h5ad',
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
 '51_B03605C2E5_WT202406020126.h5ad',
 '55_B03613E3G6_WT202403310069.h5ad',
 '59_B03612E4G6_WT202403310059.h5ad',
 '63_B03606C1E3_WT202403310061.h5ad',
 '67_A03595A1D3_WT202403310062.h5ad',
 '71_A03595A4D6_WT202403310063.h5ad',
 '75_D03468D1E3_WT202403310066.h5ad',
 '80_D03473D4E6_WT202403310070.h5ad',
 '84_B03423D1E3_WT202403310065.h5ad',
 '89_D03466A1C3_WT202403310058.h5ad',
 '99_D03470D1E3_WT202404290071.h5ad',
 '104_B03615F3G5_WT202405020035.h5ad',
 '105_D03468A4C6_WT202403310067.h5ad',
 '106_B03610F5G6_WT202403310041.h5ad',
 '108_B03423F5G6_WT202403310082.h5ad',
 '110_A03587D5E6_WT202403310038.h5ad',
 '113_A03589A5C6_WT202403310040.h5ad',
 '114_A03590A1C2_WT202403310039.h5ad',
 '115_B03609F5G6_WT202403310081.h5ad',
 '116_B03610D5E6_WT202403310078.h5ad',
 '117_B03614F1_WT202403110556.h5ad',
 '118_B03614G3_WT202403280407.h5ad',
 '119_B03616F1_WT202403280411.h5ad',
 'Y00547PC_WT202407282759.h5ad',
    
 'A03587A5C6_WT2024071215080.h5ad',
    
 'A03988A1C2_WT202407161208.h5ad',
 # 'A03591D4E5_WT2024071215074.h5ad',
 'A03590A3D6_WT202407192652.h5ad',
 'B03618D3F6_WT202407152793.h5ad',
 'B03607C4E6_WT2024071214941.h5ad',
 'A03994F1G2_WT2024071215067.h5ad',
 'A03588A1C2_WT202407161185.h5ad',
]




fig = plt.figure(figsize=(32, 12))  

gs = GridSpec(7, 8, figure=fig)  
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

plt.savefig('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/dmt_leiden_all_20250108_1.png', dpi = 450, bbox_inches = 'tight')
plt.close()


adata.write('/data/work/05.cluster/FuseMap/0106/cortex_p10_latent_embeddings_all_spatial_pretrain/dmt_leiden_20250108_1.h5ad')