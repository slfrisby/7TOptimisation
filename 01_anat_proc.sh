#!/bin/bash
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
#adds appropriate tools and options - no need to change if you have access to /imaging/mlr_imaging as the tools are in my folder 'AH', which is accessable to all users
export PATH=$PATH:/group/mlr-lab/AH/Projects/toolboxes/afni/v18.3.03
export PATH=$PATH:/imaging/local/software/anaconda/latest/x86_64/bin/
export PATH=$PATH:/group/mlr-lab/AH/Projects/toolboxes/apps/bin
export PATH=/imaging/local/software/mrtrix/v3.0.3_v2/bin/:$PATH
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
export DYLD_FALLBACK_LIBRARY_PATH=/group/mlr-lab/AH/Projects/toolboxes/afni/v18.3.03
export AFNI_3dDespike_NEW=YES
FSLDIR=/imaging/local/software/fsl/latest/x86_64/fsl
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH
export PATH=/imaging/local/software/centos7/ants/bin/ants/bin/:$PATH
export ANTSPATH=/imaging/local/software/centos7/ants/bin/ants/bin/
FSLOUTPUTTYPE=NIFTI_GZ

dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
work=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work
cd $dirp

ids=$1

#run fmriprep singularity - no freesurfer recon-all (which is for surface processing), anat-only as fMRI coreg to T1 failed multiple times
singularity run --cleanenv -B /imaging/local/software/freesurfer/7.1.1/license.txt:/opt/freesurfer/license.txt -B $dirp:/base /imaging/local/software/singularity_images/fmriprep/fmriprep-21.0.1.simg /base/data /base/derivatives/fmriprep participant --participant-label sub-"$ids" -w /base/work --fs-no-reconall --fs-license-file /opt/freesurfer/license.txt --output-spaces MNI152NLin2009cAsym:res-native --anat-only --skip-bids-validation



