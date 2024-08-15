#!/bin/bash

# makes a mean EPI of all participants' SESB to enable plotting of effect on top of EPI data. 

# add path
FSLDIR=/imaging/local/software/fsl/latest/x86_64/fsl
# this line simply configures FSL
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH

# concatenate the first 2 mean EPIs together
fslmerge -t /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work/meanEPI.nii.gz /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/halaiprep/sub-001/func/sub-001_acq-SESB_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/halaiprep/sub-002/func/sub-002_acq-SESB_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz

# add all the others on
for ids in 003 005 006 007 008 009 010 011 012 014 015 016 017 018 019; do

fslmerge -t /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work/meanEPI.nii.gz /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work/meanEPI.nii.gz /imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-SESB_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 

done

# calculate mean
fslmaths $workcond/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work/meanEPI.nii.gz -Tmean $workcond/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work/meanEPI.nii.gz
