# Warp_V-Norm

## Description

This script calculates the Warp V-Norm and the Warp V-Norm profile, as described in the paper
xxxx.
The objective is to quantify structural alterations between two analogous fiber bundles by computing diffeomorphic and rigid registration begween binar masks and mapping the resulting Warp V-Norm to the template bundle's skeleton profile.

It analyzes the structural differences between two fiber bundles assumed to be anatomically analogous. The script requires two NIfTI files representing the binary masks of these bundles (values 1 where the bundle exists, 0 otherwise) and a streamline file (.trk) representing the skeleton of the **template (fixed) bundle**. 

All files should be already aligned in the same space with an affine transformation based on anatomical images.

**Important Notes:**
*  The skeleton provided should be defined along the trajectory of the **fixed** bundle.
* All files should be already aligned in the same space with a transformation based on corresponding anatomical images.


The Warp V-Norm profile is computed through the following steps:

1.  **Registration:** Estimates the optimal rigid and diffeomorphic transformation to align the **moving** bundle mask with the **fixed** bundle mask.
2.  **Warp V-Norm Calculation:** Computes the norm-2 of the displacement vector field resulting from the total (rigid + diffeomorphic) transformation at each voxel. Higher values in the resulting Warp V-Norm map indicate greater structural divergence of the moving bundle from the fixed bundle in that spatial location.
3.  **Profile Mapping:** Maps the computed Warp V-Norm map onto the points of the provided skeleton streamline (.trk file) of the **fixed** bundle. This generates a profile of the Warp V-Norm along the trajectory of the template bundle's skeleton, providing a one-dimensional representation of the structural alterations along the bundle's extent.

This profile allows for the visualization and quantification of how the moving bundle deviates structurally from the fixed bundle along its anatomical pathway.

This metric was designed to estimate the differences between:

* A subject's fiber bundle (considered as the **moving** bundle)
* And a healthy template bundle (considered the **fixed** bundle)

## Usage

To compute Warp V Norm by computing the registration between the two bundles, and the norm-2 of the vectorial filed of the estimated transformation in each voxel:
```bash
python compute_WarpVNorm.py <bundle_fixed> <bundle_moving>  <out_dir> 
```


To compute Warp V-Norm as in previous point, but additionally this command also maps   Warp V-Norm values along the bundle-skeleton points:
```bash
python compute_WarpVNorm_profile_by_reg.py <bundle_fixed> <bundle_moving> <path_skeleton_trk> <n_points_profile> <out_dir> 
```

To compute Warp V Norm profile by using an already saved to file Warp V Norm map. In this case the registration is not computed since Warp V Norm was already estimated. The function simply maps the Warp V-Norm values along the skeleton points:
```bash
python compute_WarpVNorm_profile_by_nii.py <path_warp_vnorm_nii> <path_skeleton_trk> <n_points_profile> <path_warp_v_norm_profile_csv> 
```

