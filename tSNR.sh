#!/bin/bash
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
#adds appropriate tools and options
export PATH=$PATH:/imaging/local/software/afni/v18.3.03
export LD_LIBRARY_PATH=/imaging/local/software/afni/v18.3.03:$LD_LIBRARY_PATH

dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
work=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work
cd $dirp

# make output directory
mkdir -p $dirp/derivatives/tSNR/sub-"$ids"/

# for each run
for data in SESB SEMB SESBernst SEMBernst; do

# calculate mean and standard deviation images (for inspection). Providing no input to 3dTstat calculates the mean
3dTstat -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/"$data"_mean.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 
3dTstat -stdevNOD -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/"$data"_stdev.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 
# calculate tSNR
3dTstat -tsnr -overwrite -prefix $dirp/derivatives/tSNR/sub-"$ids"/"$data"_tSNR.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz 

done

