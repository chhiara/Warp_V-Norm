#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 17 22:49:55 2025


@author: chhiara

Script that:
     compute Warp V Norm, by computing registration RIGID and DIFF between bundles binary masks  
     then the Warp V Norm is mapped to skeleton streamlines so that Warp V Norm profile is computed

"""

import argparse
import importlib.util
import sys
import os

script_path = os.path.realpath(__file__)
script_dir = os.path.dirname(script_path)
sys.path.append(script_dir)

print(script_dir)

from WarpVNorm_library import get_WarpVnorm_Profile_computeReg

def main():
    parser = argparse.ArgumentParser(
        description='''Calculates the Warp V-Norm profile to quantify structural alterations between two analogous fiber bundles by computing diffeomorphic and rigid registration and mapping the resulting Warp V-Norm to the bundle's skeleton profile.

        This script analyzes the structural differences between two bundles
         assumed to be anatomically analogous. It requires two NIfTI files representing the binary masks of these bundles 
         (values 1 where the bundle exists, 0 otherwise) and a streamline file (.trk) representing the skeleton of the template bundle.

         NB: The skeleton should be defined along the trajectory of the fixed bundle.
        
      
        The Warp V-Norm profile is computed through the following steps:

        1.  **Registration:** Estimates the optimal rigid and diffeomorphic transformation to align the moving bundle mask  with the fixed bundle mask .
        2.  **Warp V-Norm Calculation:** Computes the norm-2 of the displacement vector field resulting from the total (rigid + diffeomorphic) transformation at each voxel. Higher values in the resulting Warp V-Norm map indicate greater structural divergence of the moving bundle from the fixed bundle in that spatial location. 
        3.  **Profile Mapping:** Maps the computed Warp V-Norm map onto the points of the provided skeleton streamline (.trk file). This generates a profile of the Warp V-Norm along the trajectory of the bundle's skeleton, providing a one-dimensional representation of the structural alterations along the bundle's extent.

        This profile allows for the visualization and quantification of how the moving bundle deviates structurally from the fixed bundle along its anatomical pathway. 
        

        This metric was designed to estimate the differences between:
        - subject's fiber bundle (considered as moving bundle) 
        - and a healthy template bundle (considered the fixed bundle)
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('bundle_fixed', help='Path to the fixed bundle file. (e.g. .nii, .nii.gz).')
    parser.add_argument('bundle_moving', help='Path to the moving bundle file. (e.g., .nii, .nii.gz).')
    parser.add_argument('path_skeleton_trk', help='Path to the skeleton bundle file. (.trk)')
    parser.add_argument('n_points_profile', help='Number of points of the skeleton profile')
    parser.add_argument('out_dir', help='Path to the output directory where the metrics will be saved.')

    args = parser.parse_args()


    # Call the function with the provided arguments
    get_WarpVnorm_Profile_computeReg(args.bundle_fixed, args.bundle_moving, args.path_skeleton_trk, args.n_points_profile, args.out_dir )
    print(f"Function 'get_WarpVnorm_Profile_computeReg' executed successfully. Check the output in: {args.out_dir}")

if __name__ == "__main__":
    main()