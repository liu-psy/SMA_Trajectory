# Spin permutation test
# t-values from APC and t-values normative models
import os
import pandas as pd
import nibabel as nib
from nilearn.datasets import fetch_atlas_surf_destrieux
from neuromaps import images, nulls
from neuromaps import parcellate
from neuromaps import stats

# set working dir
save_dir = "H:/ABCD/Release5.1/results/25_ABCD_SMA"
os.makedirs(save_dir, exist_ok=True)
os.chdir(save_dir)

# load data
APC_normative_t = pd.read_csv("APC_normative_tmaps.csv")

# extract t-maps
APC_12_t = APC_normative_t.iloc[:, 0]
APC_13_t = APC_normative_t.iloc[:, 1]
Normative_12_t = APC_normative_t.iloc[:, 2]
Normative_13_t = APC_normative_t.iloc[:, 3]

# get Destrieux atlas
destrieux = fetch_atlas_surf_destrieux()

# brain labels
labels = destrieux['labels']

parc_left = images.construct_shape_gii(destrieux['map_left'], labels=labels, intent='NIFTI_INTENT_LABEL')
parc_right = images.construct_shape_gii(destrieux['map_right'], labels=labels, intent='NIFTI_INTENT_LABEL')

# save nii
lh_path = os.path.join(save_dir, "destrieux_lh_label.gii")
rh_path = os.path.join(save_dir, "destrieux_rh_label.gii")

parc_left.to_filename(lh_path)
parc_right.to_filename(rh_path)

# relabel
parcellation = images.relabel_gifti((lh_path, rh_path), background=['Medial_wall'])
parcellation_paths = (lh_path, rh_path)

# Spatial Corelation -----------------------------------------------------------
# Persistently Low vs. Decreasing
rotated_12 = nulls.alexander_bloch(APC_12_t, atlas='fsaverage', density='10k',
                                n_perm = 10000, seed=1234, parcellation=parcellation_paths)
corr, pval = stats.compare_images(APC_12_t, Normative_12_t, nulls = rotated_12)
print(f'r = {corr:.4f}, p = {pval:.4}')

# Persistently Low vs. Increasing
rotated_13 = nulls.alexander_bloch(APC_13_t, atlas='fsaverage', density='10k',
                                n_perm = 10000, seed=1234, parcellation=parcellation_paths)
corr, pval = stats.compare_images(APC_13_t, Normative_13_t, nulls = rotated_13)
print(f'r = {corr:.4f}, p = {pval:.4}')
