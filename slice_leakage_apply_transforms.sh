#!/bin/bash
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
#adds appropriate tools and options
# no need to change if you have access to /imaging/mlr_imaging as the tools are in the folder 'AH', which is accessible to all users
FSLDIR=/imaging/local/software/fsl/latest/x86_64/fsl
# this line simply configures FSL
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
export PATH=/imaging/local/software/centos7/ants/bin/ants/bin/:$PATH
export ANTSPATH=/imaging/local/software/centos7/ants/bin/ants/bin/
FSLOUTPUTTYPE=NIFTI_GZ

#conda enviroment includes tedana toolkit and heudiconv toolkit
#need to set up if not already done/visible on your space
#source activate SLFRISBY

dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
work=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work
cd $dirp

#set ids variable (ONLY IF STEPPING THROUGH CODE)
#ids=$1

for peak in 1 2 3 4 5; do

# calculate mean MEMB functional image in native space
workcond="$work"/slice-leakage-sub-"$ids"/MEMB
fslmaths $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-MEMB_rec-t2star_run-01_space-native_desc-preproc_bold.nii.gz -Tmean $workcond/mean

# transform manually-created peak ROI from MNI to T1 space and then from T1 space to native space
antsApplyTransforms -d 3 -i $dirp/scripts/slice_leakage/peak$peak.nii.gz -r $workcond/mean.nii.gz -o $work/slice-leakage-sub-"$ids"/MEMB/peak$peak.nii.gz --default-value 0 --float 1 -n NearestNeighbor --transform [$dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-01_from-T1w_to-MEMB_mode-image_xfm.mat,0] --transform $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5 --transform identity --transform identity

done



