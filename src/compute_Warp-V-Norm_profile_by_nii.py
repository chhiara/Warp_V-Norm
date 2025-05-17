#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 17 22:49:55 2025


@author: chhiara

Script that:
     compute Warp V Norm profile, by using an already computed Warp V Norm map.
      The registration RIGID and DIFF between bundles binary masks  
     was already performed to compute the Warp V Norm, and is not performed by this script.
     This script simply map the  Warp V Norm values to the skeleton streamlines to compute the Warp V Norm profile

"""

import argparse
import importlib.util
import sys
import os
from Warp_V-Norm_library import get_WarpVnorm_Profile_from_warpVnorm

def main():
    parser = argparse.ArgumentParser(
        description='''Calculates the Warp V-Norm profile to quantify structural alterations between two analogous fiber bundles by using an already saved Warp V Norm map, that is mapped the bundle's skeleton profile.

        This script analyzes the structural differences between two bundles
         assumed to be anatomically analogous. It requires the Warp V-Norm, already saved to file 
         and a streamline file (.trk) representing the skeleton of the template bundle.

         NB: The skeleton should be defined along the trajectory of the fixed bundle, used in the registration to compute Warp V-Norm
        
        The Warp V-Norm profile is computed by mapping the Warp V-Norm onto the points of the provided skeleton streamline (.trk file). This generates a profile of the Warp V-Norm along the trajectory of the bundle's skeleton, providing a one-dimensional representation of the structural alterations along the bundle's extent.

        This profile allows for the visualization and quantification of how the moving bundle deviates structurally from the fixed bundle along its anatomical pathway. 
        

        This metric was designed to estimate the differences between:
        - subject's fiber bundle (considered as moving bundle) 
        - and a healthy template bundle (considered the fixed bundle)
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('path_warp_vnorm_nii', help='Path to the mWarp V-Norm. (e.g., .nii, .nii.gz).')
    parser.add_argument('path_skeleton_trk', help='Path to the skeleton bundle file. (.trk)')
    parser.add_argument('n_points_profile', help='Number of points of the skeleton profile')
    parser.add_argument('path_warp_v_norm_profile_csv', help='Path to the csv file that where the Warp V-Norm profile will be saved.')

    args = parser.parse_args()


    # Check if all compulsory arguments are provided (argparse usually does this,
    # but let's be explicit for clarity in this context)
    if not all([args.bundle_fixed, args.bundle_moving, args.path_skeleton_trk, args.n_points_profile, args.out_dir]):
        parser.print_help()
        sys.exit(1)

    # Call the function with the provided arguments
    get_WarpVnorm_Profile_from_warpVnorm(args.path_skeleton_trk, args.path_warp_vnorm_nii, args.path_warp_v_norm_profile_csv, args.n_points_profile)
    print(f"Function 'get_WarpVnorm_Profile_computeReg' executed successfully. Check the output in: {args.out_dir}")

if __name__ == "__main__":
    main()