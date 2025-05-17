#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 13 15:49:55 2025


@author: chhiara

Script to run RIGID - DIFF registration between bundles binary masks to compute Warp V Norm metrics 

"""
#%%
import nibabel as nib
import numpy as np
import os
import errno
import sys
from pathlib import Path
import pandas as pd
import dipy.stats.analysis as dsa
import dipy.tracking.streamline as dts
from dipy.io.streamline import load_trk

#%%
def print_run_bash(bash_cmd):
    print(bash_cmd)
    os.system(bash_cmd)
    print()

def dice_score(y_gt,y_pred):
    """   
    Parameters
    ----------
    y_gt : np.array or iterable.
        vector of the ground thruth.
    y_pred : np.array or iterable.
        vectors of the binary prediction.

    Returns
    -------
    The Dice Score Coefficient.

    """
    assert len(y_pred) == len(y_gt), "y_pred: {len(y_pred)}    y_gt: {len(y_gt)}. They should have instead the same size! "
    y_pred= y_pred if type(y_pred) == np.array else np.array(y_pred)
    y_gt= y_gt if type(y_gt) == np.array else np.array(y_gt)
    
    #compute Dice score
    #1)compute the intersection between ground truth and prediction

    y_pred = (np.rint(y_pred)).astype(int)
    y_gt = (np.rint(y_gt)).astype(int)

    pos_pred_test=(y_pred==1)
    pos_gt_test=(y_gt==1)

    intersec=np.sum((pos_pred_test)*(pos_gt_test))
    
    #compute the noramlization factor
    norm_factor=np.sum(y_pred)+np.sum(y_gt)
    #dice score
    d_score=(2*intersec)/norm_factor
    return(d_score)


def symlink_force(target, link_name):
    try:
        os.symlink(target, link_name)
    except OSError as e:
        if e.errno == errno.EEXIST:
            os.remove(link_name)
            os.symlink(target, link_name)
        else:
            raise e


def save_info_reg(bundle_fixed, bundle_moving, ants_call, output_prefix_reg ):
    dir_reg= os.path.dirname(output_prefix_reg)
    with open(f"{dir_reg}/log_reg.txt", "w") as log_f:
        log_f.write(f"bundle_fixed  {bundle_fixed}\n\n")
        log_f.write(f"bundle_moving  {bundle_moving}\n\n")
        log_f.write(ants_call)
    symlink_force(bundle_fixed, (dir_reg + "/" + "fixed_" + os.path.basename(bundle_fixed)) )
    symlink_force(bundle_moving, (dir_reg + "/" +  "moving_" + os.path.basename(bundle_moving)) )


def rigid_reg_bdls(bundle_fixed, bundle_moving, output_prefix_reg):
        
    #check if the images are in the same space
    m = nib.load(bundle_moving)
    f = nib.load(bundle_fixed)
    
    assert (f.affine == m.affine).all() , f"Errore the affine of {bundle_fixed} and {bundle_moving} are not the same. Register the two images to be in same space!"
        
    #ants call analogous to antsRegistrationSyN.sh -t a , but that uses metrics_MSQ and not MI, and only the rigid part is selected
    ants_call=f"antsRegistration --verbose 1 -d 3 --float 0 --collapse-output-transforms 1 -o [{output_prefix_reg},{output_prefix_reg}Warped.nii.gz,{output_prefix_reg}InverseWarped.nii.gz] -n NearestNeighbor --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] --initial-moving-transform [{bundle_fixed},{bundle_moving},1] -t Rigid[0.1] -m MeanSquares[{bundle_fixed},{bundle_moving},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox"
    
    save_info_reg(bundle_fixed, bundle_moving, ants_call, output_prefix_reg )

    path_bdl_reg=f"{output_prefix_reg}Warped.nii.gz"
    if not os.path.exists(path_bdl_reg):
        print_run_bash(ants_call)
    else:
        print(f"{path_bdl_reg} already exists. Registration Skipped")

def affine_reg_bdls(bundle_fixed, bundle_moving, output_prefix_reg):
        
    #check if the images are in the same space
    m = nib.load(bundle_moving)
    f = nib.load(bundle_fixed)
    
    assert (f.affine == m.affine).all() , f"Errore the affine of {bundle_fixed} and {bundle_moving} are not the same. Register the two images to be in same space!"
    
        
    #bash_command=f"antsRegistrationSyN.sh -d 3 -f {bundle_fixed} -m {bundle_moving} -o {output_prefix_reg} -n 1 -t a"
    
    #ants call analogous to antsRegistrationSyN.sh -t a , but that uses metrics_MSQ and not MI
    ants_call=f"antsRegistration --verbose 1 -d 3 --float 0 --collapse-output-transforms 1 -o [{output_prefix_reg},{output_prefix_reg}Warped.nii.gz,{output_prefix_reg}InverseWarped.nii.gz] -n NearestNeighbor --use-histogram-matching 0 --winsorize-image-intensities [0.005,0.995] --initial-moving-transform [{bundle_fixed},{bundle_moving},1] -t Rigid[0.1] -m MeanSquares[{bundle_fixed},{bundle_moving},1,32,Regular,0.25] --convergence [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox -t Affine[0.1] -m MeanSquares[{bundle_fixed},{bundle_moving},1,32,Regular,0.25] -c [1000x500x250x100,1e-6,10] --shrink-factors 8x4x2x1 --smoothing-sigmas 3x2x1x0vox"
    
    save_info_reg(bundle_fixed, bundle_moving, ants_call, output_prefix_reg )

    path_bdl_reg=f"{output_prefix_reg}Warped.nii.gz"
    if not os.path.exists(path_bdl_reg):
        print_run_bash(ants_call)
    else:
        print(f"{path_bdl_reg} already exists. Registration Skipped")

def nl_reg_bdls(bundle_fixed, bundle_moving, bundle_moving_affined, mat_affine, output_prefix_reg):
        
    #check if the images are in the same space
    m = nib.load(bundle_moving_affined)
    f = nib.load(bundle_fixed)
    
    assert (f.affine == m.affine).all() , f"Errore the affine of {bundle_fixed} and {bundle_moving} are not the same. Register the two images to be in same space!"
    
    grad_step=0.25  #gradient step of the 1',2' and 3' pyramidal levels
    its1="5000x0x0"     #iterations of the 1' pyramidal level
    its23="0x100x100"    #iterations of the 2', 3' pyramidal levels

    #compute nl registration using as initial image bdl already affined
    ants_call = f"antsRegistration -d 3 --float 1  -v 1  -m MeanSquares[{bundle_fixed},{bundle_moving_affined},1] -c [{its1},-0.01,5 ] -t SyN[{grad_step}] -s 1x0.5x0vox -f 4x2x1 -l 1 -z 1  -m MeanSquares[{bundle_fixed},{bundle_moving_affined},1] -c [{its23},-0.01,5 ] -t SyN[{grad_step}] -s 1x0.5x0vox -f 4x2x1 -l 1 -z 1  -o [{output_prefix_reg}, {output_prefix_reg}Warped.nii.gz, {output_prefix_reg}InverseWarped.nii.gz] -n NearestNeighbor"
    
    #re-compute warped bundle applying sum of affine and warp-nl. In this way the resulting warped bdl is obtained with only 1 interpolation.
    #  THis result will overwrite the one directly obtained from registration
    warp_nl=f"{output_prefix_reg}0Warp.nii.gz"
    bdl_warped_nl=f"{output_prefix_reg}Warped.nii.gz"

    warp_bdl_affinewarp_1interpolation=f"antsApplyTransforms -d 3 -i {bundle_moving} -r {bundle_fixed} -t {warp_nl} {mat_affine} -n NearestNeighbor -o {bdl_warped_nl} -v"
    
    #save log file of registration

    save_info_reg(bundle_fixed, bundle_moving, ants_call, output_prefix_reg )

    if not os.path.exists(bdl_warped_nl):
        print_run_bash(ants_call)
        print_run_bash(warp_bdl_affinewarp_1interpolation)
    else:
        print(f"{bdl_warped_nl} already exists. Registration Skipped")
    
def compute_dscs_before_after_reg(template_fixed, bundle_moving, bundle_moving_warped, cutoff_bin_template=0.5):
    template_np=nib.load(template_fixed).get_fdata()
    template_bin_np=template_np>cutoff_bin_template   #binarize the template to compute the DSC

    bdl_apss_np=nib.load(bundle_moving).get_fdata()
    bdl_apss_np_reg=nib.load(bundle_moving_warped).get_fdata()

    dsc_before_reg=round(dice_score(template_bin_np, bdl_apss_np),3)
    dsc_after_reg=round(dice_score(template_bin_np, bdl_apss_np_reg),3)
        
    print(f"DSC before reg: {dsc_before_reg}")
    print(f"DSC after reg: {dsc_after_reg}")
    return dsc_before_reg, dsc_after_reg

def get_info_from_affine(affine_ants_mat_path):
    #1. transform mat -> txt
    affine_ants_txt_path=affine_ants_mat_path.replace(".mat",".txt")
    bash_cmd=f"ConvertTransformFile 3 {affine_ants_mat_path}  {affine_ants_txt_path}"
    print_run_bash(bash_cmd)

    #2. transform matrix in txt from itk format -> to mrtrix format
    affine_txt_mrtrix_path = affine_ants_txt_path.replace(".txt","_mrtrix.txt")
    bash_cmd=f"transformconvert {affine_ants_txt_path} itk_import {affine_txt_mrtrix_path} -f"
    print_run_bash(bash_cmd)

    #3. compute metrics with transformcalc
    affine_metrics_path=affine_txt_mrtrix_path.replace(".txt","_metrics.txt")
    bash_cmd=f"transformcalc {affine_txt_mrtrix_path} decompose {affine_metrics_path} -f"
    print_run_bash(bash_cmd)

    #4. from the file obtained with transformcalc get the metrics
    with open(affine_metrics_path, "r") as metrics_file:
        lines_metrics_file=metrics_file.readlines()
    #print(lines_metrics_file)
    jac_det=round(float((([l for l in lines_metrics_file if "jacobian_det:" in l][0]).replace("jacobian_det: ", "")).replace("\n","")),3)

    vect_translation=((([l for l in lines_metrics_file if "translation" in l][0]).replace("translation: ", "")).replace("\n","")).split(" ")
    vect_translation_np=np.array([round(float(n),3) for n in vect_translation ])
    norm_transl=round(np.linalg.norm(vect_translation_np),3)


    return norm_transl, jac_det



def create_warp_affine(bundle_fixed, mat_affine, warp, warp_affine):
    if os.path.exists(warp_affine):
        print(f"{warp_affine} exists. Creation of the total warp skipped")
    else:
        #combine affine and warp in a unique warp_field saved to file as warp_affine
        ants_call=f"antsApplyTransforms -d 3 -t {warp} {mat_affine} -o [{warp_affine},1] -r {bundle_fixed} -v"
        print_run_bash(ants_call)



def compute_displacement_tot(bundle_fixed, bundle_moving_warped, mat_affine, warp, warp_affine, warp_v_norm_img):

    #combine affine and warp in a unique warp_field saved to file as warp_affine
    create_warp_affine(bundle_fixed, mat_affine, warp, warp_affine)

    #compute sum(norm(vector_filed)) over all the warped coordinates of the warped-bundle in the fixed space
    warp_affine_np=nib.load(warp_affine).get_fdata()
    bundle_moving_warped_img = nib.load(bundle_moving_warped)
    bundle_moving_warped_np = bundle_moving_warped_img.get_fdata()
    assert np.sum(bundle_moving_warped_np>0) == np.sum(bundle_moving_warped_np==1)

    #select the points of the warp where the warped-bundle is defined (voxel ==1 and not zero)
    #this are the interesting points for which I compute the vectorial field
    warped_bdl_true=np.nonzero(bundle_moving_warped_np)
    
    norm_vect_map=np.zeros(bundle_moving_warped_np.shape)
    norms=[]

    n_true_warped_vox=len(warped_bdl_true[0])
    i=0
    for i in range(n_true_warped_vox):
        x = warped_bdl_true[0][i]
        y = warped_bdl_true[1][i]
        z = warped_bdl_true[2][i]
        
        v=np.array([warp_affine_np[x,y,z,0,0], #x componenet of the dispalcement-vectors in the fixed space 
                    warp_affine_np[x,y,z,0,1], #y componenet of the dispalcement-vector
                    warp_affine_np[x,y,z,0,2]]) #z componenet of the dispalcement-vector


        normm = np.linalg.norm(v)
        norms.append(normm)
        norm_vect_map[x,y,z]=normm

    out_nib=nib.Nifti1Image(norm_vect_map, affine=bundle_moving_warped_img.affine, header=bundle_moving_warped_img.header)
    out_nib.to_filename(warp_v_norm_img)
    print(f"{warp_v_norm_img} Saved to file")   #--> this will store the local measure of displacement

    mean_warp_v_norm_img = round(np.mean(np.array(norms)),3) #mean of the norm of the displacement vectors--> this will be an aggregate measure of displacement
    
    return mean_warp_v_norm_img

def compute_JacDet_tot(bundle_fixed, bundle_moving_warped, mat_affine, warp, warp_affine, jac_det_warp_affine):

    #create combined transformation, combining affine and warp if it does not exists
    create_warp_affine(bundle_fixed, mat_affine, warp, warp_affine)

    #create determinant of the jacobian
    ants_call=f"CreateJacobianDeterminantImage 3 {warp_affine} {jac_det_warp_affine}"
    print_run_bash(ants_call)

    #compute the mean of the jacobian in the coordinates of the warped bundles, in fixed space
    jac_det_np = nib.load(jac_det_warp_affine).get_fdata()
    bundle_moving_warped_np = nib.load(bundle_moving_warped).get_fdata()

    assert jac_det_np.shape == bundle_moving_warped_np.shape
    assert np.sum(bundle_moving_warped_np==1) == np.sum(bundle_moving_warped_np>0)

    mean_jac_det_bdl_warped = round(np.mean(jac_det_np[bundle_moving_warped_np==1]),3)

    return mean_jac_det_bdl_warped




def get_warp_vnorm_metrics_rigid_diff_reg(bundle_fixed, bundle_moving, out_dir):
    """
        Function to get metrics compute to evaluate how much a bundle is different from a reference bundle
        
        The function computes: 
            0) RIGID reg + diff registration
        
       of a bundle_moving (of a subject) --> assuming this image is binary, eg in the voxels nonly 0 or 1 are present
            versus 
        a bundle fixed (a reference, eg a template or a prob map of the bundle)  --> assuming this image is either binary, or with values in the range [0,1]
        
        The registration is useful to evaluate how much a bundle is different from a reference bundle. This is done computing:
            1)the  whole vectorial field  of the total affine + nl  transformation
    
        From this compute:
            
        2) the norm of the vectors in the area of the warped bundle (in fixed space)
                 >  average the norms'values to have a total measure of displacement -> gloabal measure
                 > save a nifti map storing the norm of the vector of displacement for each voxel of the warped bundle  --> local meaasure

        3) the Jacobian determinant of the warped bundle (in fixed space)
                > average the jacobian det in the area of the warped bundle to have a total measure of variation of volume
                     NB: if Jacdet<1  -> the template decrease the volume to match moving_bundle. So moving_bundle is expanding (and viceversa)
                 > save a nifti map storing the jacobian determinant of the vectorial field of the trasformation for each voxel of the warped bundle --> local meaasure

    Parameters
    ----------
    bundle_fixed : string
        path to the nifti image of the fixed bundle (eg the template bundle or the probability map of the bundle).
    bundle_moving : string
        path to the nifti image of the moving bundle (eg the bundle of a specific subject I want to characterize).
    out_dir : string
        path to output_directory.

    Returns
    -------
    dsc_before_aff_reg : float
        DSC of the bundle_moving versus the bundle fixed (taken as ground truth) before all the registration.
    dsc_after_aff_reg : float
        DSC of the bundle_moving (warped with affine registration) with respect to the bundle fixed (taken as ground truth).
    dsc_after_nl_reg : float
        DSC of the bundle_moving (warped with affine and diff registration) with respect to the bundle fixed (taken as ground truth).
    mean_warp_v_norm_img : float
        mean of the norm2 of the vectors (representing the whole trasformation affine + diff) in te area of the warped bundle.
    mean_jac_det_bdl_warped : float
       mean of the Determinant of the jacobian of the vectors (representing the whole trasformation affine + diff) in the area of the warped bundle.

    """
    
    b_moving_name=os.path.basename(bundle_moving).split(".")[0]
    
    #define and create out directories
    dir_reg_aff=f"{out_dir}/rigid/"
    dir_reg_nl=f"{out_dir}/non-linear/"
    dir_metrics_reg_nii=f"{out_dir}/reg-metrics/"
    
    for dir_path in [dir_reg_aff, dir_reg_nl, dir_metrics_reg_nii]:
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
    
    #path rigid reg
    output_prefix_reg_aff= f"{dir_reg_aff}/{b_moving_name}_rig_"
    bundle_moving_affined=f"{output_prefix_reg_aff}Warped.nii.gz"
    mat_affine=f"{output_prefix_reg_aff}0GenericAffine.mat"
    
    #path non linear reg
    output_prefix_reg_nl= f"{dir_reg_nl}/{b_moving_name}_diff_"
    bundle_moving_nl_warped=f"{output_prefix_reg_nl}Warped.nii.gz"
    warp=f"{output_prefix_reg_nl}0Warp.nii.gz"
    
    #path of metrics of total regisration rigid + dffeomorphic 
    warp_rigid=f"{dir_metrics_reg_nii}/rigid_diff_warp.nii.gz"
    
    warp_v_norm_img=f"{dir_metrics_reg_nii}/warp_v-norm_of_displacement_rigid-diff.nii.gz"
    jac_det_warp_rigid=f"{dir_metrics_reg_nii}/jac_det_of_displacement_rigid-diff.nii.gz"

    
    #compute rigid reg
    rigid_reg_bdls(bundle_fixed, bundle_moving, output_prefix_reg_aff)
    
    #compute nl reg
    nl_reg_bdls(bundle_fixed, bundle_moving, bundle_moving_affined, mat_affine, output_prefix_reg_nl)
    
    #compute dsc before after rigid reg 
    dsc_before_rig_reg, dsc_after_rig_reg = compute_dscs_before_after_reg(bundle_fixed, bundle_moving, bundle_moving_affined)
    
    #compute dsc before after nl reg
    dsc_before_nl_reg, dsc_after_nl_reg = compute_dscs_before_after_reg(bundle_fixed, bundle_moving_affined, bundle_moving_nl_warped)
    
    assert dsc_before_nl_reg == dsc_after_rig_reg
    
    #compute vector displacement of the registration rigid + diff --> compute the norm of the vector displacement
    mean_warp_v_norm_img = compute_displacement_tot(bundle_fixed, bundle_moving_nl_warped, mat_affine, warp, warp_rigid, warp_v_norm_img)
    
    #compute jacobian tot
    mean_jac_det_bdl_warped = compute_JacDet_tot(bundle_fixed, bundle_moving_nl_warped, mat_affine, warp, warp_rigid, jac_det_warp_rigid)
    
    return dsc_before_rig_reg, dsc_after_rig_reg, dsc_after_nl_reg, mean_warp_v_norm_img, mean_jac_det_bdl_warped



def get_WarpVnorm_Profile_from_warpVnorm(path_skeleton_trk, path_warp_vnorm_nii, path_warp_v_norm_profile_csv, n_points_profile):

        skeleton=(load_trk(path_skeleton_trk, "same").streamlines)
        
        #get the warp vnorm map
        warp_v_norm_img = nib.load(path_warp_vnorm_nii)
        warp_v_norm_img_np = warp_v_norm_img.get_fdata()
        
     
        reference_streamline = skeleton[0]
        
        oriented_skeleton = dts.orient_by_streamline(skeleton, reference_streamline)
        #-------------------------------

        #compute skeleton profile of vnorm
        skeleton_profile = dsa.afq_profile(warp_v_norm_img_np, oriented_skeleton,  warp_v_norm_img.affine, n_points=n_points_profile)
        skeleton_profile = [round(x,3) for x in skeleton_profile ]
        
        
        #save the complete data per subjects and for all subjects together
        data_dict_complete={
                    "skeleton_v_norm":skeleton_profile,
                    "point": [ i for i in range(0,n_points_profile)]}
        
        skeleton_vnorm = pd.DataFrame(data_dict_complete)    


    
        #save the whole profile in dataframes converted in csvs
        skeleton_vnorm.to_csv(path_warp_v_norm_profile_csv, index=False)
     
     
def get_WarpVnorm_Profile_computeReg(bundle_fixed, bundle_moving, path_skeleton_trk, n_points_profile, out_dir ):
    
    get_warp_vnorm_metrics_rigid_diff_reg(bundle_fixed, bundle_moving, out_dir)
    
    dir_metrics_reg_nii=f"{out_dir}/reg-metrics/"
        
    path_warp_vnorm_nii=f"{dir_metrics_reg_nii}/warp_v-norm_of_displacement_rigid-diff.nii.gz"
    path_warp_v_norm_profile_csv = f"{dir_metrics_reg_nii}/warp_v-norm_of_displacement_rigid-diff__profile_skeleton.trk"
    
    get_WarpVnorm_Profile_from_warpVnorm(path_skeleton_trk, path_warp_vnorm_nii, path_warp_v_norm_profile_csv, n_points_profile)