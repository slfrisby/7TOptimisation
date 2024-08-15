#!/bin/bash
echo "
++++++++++++++++++++++++" 
echo +* "Set up script run environment" 
#adds appropriate tools and options - no need to change if you have access to /imaging/mlr_imaging as the tools are in the folder 'AH', which is accessible to all users
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

#conda enviroment includes tedana toolkit and heudiconv toolkit
#need to set up if not already done/visible on your space
source activate SLFRISBY

dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
cd $dirp

$dirp/scripts/./sub_7Tpilot_ME_slice_leakage.sh "$ids"

$dirp/scripts/./sub_7Tpilot_SE_slice_leakage.sh "$ids"

