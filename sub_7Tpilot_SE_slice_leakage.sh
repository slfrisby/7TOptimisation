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
# this line simply configures FSL
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
work=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work
cd $dirp

#set ids variable; comment out if using sub_job.sh script to initiate this script
#ids=$1

#####DO MP2RAGE PROCESSING#####
if [ ! -f $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_desc-preproc_T1w.nii.gz ]; then
echo "MP2RAGE processed files missing - running preprocessing steps"
#call script that creates denoised MP2RAGE and run CAT12 segment. Comment out for sub001 because of missing reverse phase file for ptx
$dirp/scripts/./convert_to_BIDS.sh "$ids"
#run fmriprep anat-only as func runs coreg fails randomly
$dirp/scripts/./01_anat_proc.sh "$ids"
else
echo "MP2RAGE processed files found - skipping preprocessing"
fi


#####RUN fMRI processing#####
for data in SESB SEMB; do

run=$dirp/data/sub-"$ids"/func/sub-"$ids"_task-semantic_acq-"$data"_run-01_bold.nii.gz

filename=$(basename $run .nii.gz)

workcond="$work"/slice-leakage-sub-"$ids"/"$data"
rm -rf $workcond/*
mkdir $workcond -p

# don't despike EPIs - this command gets rid of background noise, which we want to keep for slice leakage testing. Instead, to keep the naming conventions consistent for the rest of the code, just copy the right files in and call them the right names
cp $run $workcond/tmpdata1.nii.gz

#Slice time correction
TR=$(cat "$dirp"/data/sub-"$ids"/func/$filename.json | jq '.RepetitionTime')
# get the times of each slice, store each time on a separate line, and save as a .txt file
cat "$dirp"/data/sub-"$ids"/func/$filename.json | python2.7 "$dirp"/scripts/slice_out.py > $workcond/timings.txt
# conduct the actual slice timing correction. By default, 3dTshift searches for slice timing data in the header of the file it is trying to correct. Our scans have that information in the .json file instead, which is we provide a separate text file.
3dTshift -ignore 0 -prefix $workcond/atmpdata1.nii.gz -tpattern @$workcond/timings.txt -TR $TR $workcond/tmpdata1.nii.gz

#Realign EPIs using 1st echo and apply to other echoes
# this is a calculation, but all it does is get the first volume, do nothing to it (-expr 'a') and save it 
3dcalc -a $workcond/atmpdata1.nii.gz[0]  -expr 'a' -prefix $workcond/tmpbase.nii.gz 
# align subsequent volumes to the first volume 
3dvolreg -overwrite -prefix $workcond/ratmpdata1.nii.gz -base $workcond/tmpbase.nii.gz -dfile $workcond/ratmpdata.1D -1Dmatrix_save $workcond/ratmpdata.aff12.1D $workcond/atmpdata1.nii.gz
# get the 6 rigid-body motion parameters from the .1D file created above and save those to a new .1D file
1dcat $workcond/ratmpdata.1D[1..6] > $workcond/"$filename"_motion.1D 
# save that same information as a text file
cat $workcond/"$filename"_motion.1D > $workcond/"$filename"_motion.txt
# do the final realignment using an affine transformation (this is what 3dAllineate does) and using the 6 rigid-body motion parameters generated earlier
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply $workcond/ratmpdata.aff12.1D  -base $workcond/tmpbase.nii.gz -input $workcond/atmpdata1.nii.gz -prefix $workcond/ratmpdata1.nii.gz -overwrite

#Setup files for TOPUP
# construct merged file contains 5 volumes of AP phase encoding direction and 5 volumes of PA
fslmerge -t $workcond/merged.nii.gz $dirp/data/sub-$ids/fmap/*"$data"*AP*.nii.gz $dirp/data/sub-$ids/fmap/*"$data"*PA*.nii.gz
# setup acquisition parameter file. Get total readout time (time from centre of first echo to centre of last, in seconds)
totalreadouttime=$(cat "$dirp"/data/sub-"$ids"/fmap/*"$data"*AP*.json | jq '.TotalReadoutTime')
# setup first 5 lines of acquisition parameter file. It doesn't matter which of AP and PA is labelled with 1 and -1 - all that matters is that they are orthogonal.
x=(0,-1,0,$totalreadouttime)
# copy this line 5 times and append to the file
echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt
# repeat for the next 5 lines
x=(0,1,0,$totalreadouttime)
echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt;echo $x | tr ',' ' ' >> $workcond/acqp.txt

#run TOPUP...
topup --imain=$workcond/merged.nii.gz --datain=$workcond/acqp.txt --config=b02b0.cnf --out=$workcond/topup --iout=$workcond/distcorr.nii.gz

#...and apply the field map to T2star and ICA-cleaned EPIs. --inindex means that the first line of acqp.text describes the main file that we want to transform (i.e. our data is acquired with AP phase encoding).
applytopup --imain=$workcond/ratmpdata1.nii.gz --inindex=1 --datain=$workcond/acqp.txt --topup=$workcond/topup --method=jac --out=$workcond/uratmpdata1.nii.gz

#Don't coregister or normalise - if you want to see the results on a brain you should back-project the structural image into functional space. Instead, just copy the latest files over
cp $workcond/uratmpdata1.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_run-01_space-native_desc-preproc_bold.nii.gz

done




