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
#call script that creates denoised MP2RAGE and run CAT12 segment
$dirp/scripts/./00_convert_to_BIDS.sh "$ids"
#run fmriprep anat-only as func runs coreg fails randomly
$dirp/scripts/./01_anat_proc.sh "$ids"
else
echo "MP2RAGE processed files found - skipping preprocessing"
fi


#####RUN fMRI processing#####
for data in MESB MEMB; do

run=$dirp/data/sub-"$ids"/func/sub-"$ids"_task-semantic_acq-"$data"_run-01_echo-1_bold.nii.gz

filename=$(basename $run .nii.gz)

workcond="$work"/sub-"$ids"/"$data"
rm -rf $workcond/*
mkdir $workcond

#Despike EPIs
3dDespike -overwrite -prefix $workcond/tmpdata1.nii.gz $run
3dDespike -overwrite -prefix $workcond/tmpdata2.nii.gz $(echo "${run/-1_/-2_}")
3dDespike -overwrite -prefix $workcond/tmpdata3.nii.gz $(echo "${run/-1_/-3_}")

#Slice time correction
TR=$(cat "$dirp"/data/sub-"$ids"/func/$filename.json | jq '.RepetitionTime')
# get the header and save it as a text file
fslhd $run > $workcond/tmp.txt
# get the 8th line down (which gives dim3 = the dimension in the z direction = the number of slices) and save it as another text file
sed '8q;d' $workcond/tmp.txt > $workcond/tmp1.txt
# delete the first 15 characters from the line (leaving only the number of slices)
sed 's/...............//' $workcond/tmp1.txt > $workcond/tmp2.txt
# get the number of slices as a variable
slices=$(cat $workcond/tmp2.txt)
# from the .json file for the first echo scan, get the times of each slice, store each time on a separate line, and save as a .txt file. The -oe flag gets rid of the words "SliceTiming"
grep "SliceTiming" --after-context=$slices "$dirp"/data/sub-"$ids"/func/$filename.json | grep -oe '\([0-9.]*\)' | cat > $workcond/timing.txt
# conduct the actual slice timing correction. By default, 3dTshift searches for slice timing data in the header of the file it is trying to correct. Our scans have that information in the .json file instead, which is we provide a separate text file.
3dTshift -ignore 0 -prefix $workcond/atmpdata1.nii.gz -tpattern @$workcond/timing.txt -TR $TR $workcond/tmpdata1.nii.gz
3dTshift -ignore 0 -prefix $workcond/atmpdata2.nii.gz -tpattern @$workcond/timing.txt -TR $TR $workcond/tmpdata2.nii.gz
3dTshift -ignore 0 -prefix $workcond/atmpdata3.nii.gz -tpattern @$workcond/timing.txt -TR $TR $workcond/tmpdata3.nii.gz

#Realign EPIs (do 1st echo first, then apply to other echoes)
# this is a calculation, but all it does is get the first volume of the first echo, do nothing to it (-expr 'a') and save it 
3dcalc -a $workcond/atmpdata1.nii.gz[0]  -expr 'a' -prefix $workcond/tmpbase.nii.gz 
# align subsequent volumes of the first echo dataset to the first volume 
3dvolreg -overwrite -prefix $workcond/ratmpdata1.nii.gz -base $workcond/tmpbase.nii.gz -dfile $workcond/ratmpdata.1D -1Dmatrix_save $workcond/ratmpdata.aff12.1D $workcond/atmpdata1.nii.gz
# get the 6 rigid-body motion parameters from the .1D file created above and save those to a new .1D file
1dcat $workcond/ratmpdata.1D[1..6] > $workcond/"$filename"_motion.1D 
# save that same information as a text file
cat $workcond/"$filename"_motion.1D > $workcond/"$filename"_motion.txt
# do the final realignment using an affine transformation (this is what 3dAllineate does) and using the 6 rigid-body motion parameters generated earlier
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply $workcond/ratmpdata.aff12.1D  -base $workcond/tmpbase.nii.gz -input $workcond/atmpdata1.nii.gz -prefix $workcond/ratmpdata1.nii.gz -overwrite
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply $workcond/ratmpdata.aff12.1D  -base $workcond/tmpbase.nii.gz -input $workcond/atmpdata2.nii.gz -prefix $workcond/ratmpdata2.nii.gz -overwrite
3dAllineate -final wsinc5 -cubic -float -1Dmatrix_apply $workcond/ratmpdata.aff12.1D  -base $workcond/tmpbase.nii.gz -input $workcond/atmpdata3.nii.gz -prefix $workcond/ratmpdata3.nii.gz -overwrite

#t2star and tedana combination 
# get the value of each echo from the .json file and convert it from seconds to miliseconds
e1=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$filename.json | jq '.EchoTime')*1000 | bc)
tmp=$(echo "${filename/-1_/-2_}")
e2=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$tmp.json | jq '.EchoTime')*1000 | bc)
tmp=$(echo "${filename/-1_/-3_}")
e3=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$tmp.json | jq '.EchoTime')*1000 | bc)
# run brain extraction. One of the outputs is a brain mask at native EPI resolution in EPI space
bet $workcond/ratmpdata1.nii.gz $workcond/brain -m -f 0.2
# tedana workflow
tedana -d $workcond/ratmpdata1.nii.gz $workcond/ratmpdata2.nii.gz $workcond/ratmpdata3.nii.gz -e $e1 $e2 $e3 --mask $workcond/brain_mask.nii.gz --out-dir $workcond/tedana --fittype curvefit --n-threads 16 --maxit 500 --maxrestart 50 
# this just makes sure that the outputs have headers in FSL's standard space. This means that transforms will be applied correctly.
fslreorient2std $workcond/tedana/desc-optcom_bold $workcond/tedana/desc-optcom_bold
fslreorient2std $workcond/tedana/desc-optcomDenoised_bold $workcond/tedana/desc-optcomDenoised_bold

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

#transform native merged files into tedana output space (they are currently in native EPI space)
#select first volume of the merged reverse phase file
fslselectvols -i $workcond/merged -o $workcond/tedana/epi_native_vol --vols=0
#register that volume to tedana output
flirt -in $workcond/tedana/epi_native_vol -ref $workcond/tedana/desc-optcom_bold.nii.gz -out $workcond/tedana/cepi_native_vol -omat $workcond/tedana/cepi_native_vol.mat -dof 6
#apply the transform to merged files to match tedana outputs 
flirt -in $workcond/merged -ref $workcond/tedana/cepi_native_vol -applyxfm -init $workcond/tedana/cepi_native_vol.mat -interp nearestneighbour -out $workcond/tedana/merged

#run TOPUP...
topup --imain=$workcond/tedana/merged.nii.gz --datain=$workcond/acqp.txt --config=b02b0.cnf --out=$workcond/tedana/topup --iout=$workcond/tedana/distcorr.nii.gz

#...and apply the field map to T2star and ICA-cleaned EPIs. --inindex means that the first line of acqp.text describes the main file that we want to transform (i.e. our data is acquired with AP phase encoding).
applytopup --imain=$workcond/tedana/desc-optcomDenoised_bold.nii.gz --inindex=1 --datain=$workcond/acqp.txt --topup=$workcond/tedana/topup --method=jac --out=$workcond/udesc-optcomDenoised_bold.nii.gz
applytopup --imain=$workcond/tedana/desc-optcom_bold.nii.gz --inindex=1 --datain=$workcond/acqp.txt --topup=$workcond/tedana/topup --method=jac --out=$workcond/udesc-optcom_bold.nii.gz

#Coregistration. Coregister mean EPI to structural image, then apply transforms to all EPIs - keep at native resolution in MNI space (i.e., 2.5mm)
# calculate mean EPI across time
fslmaths $workcond/udesc-optcom_bold.nii.gz -Tmean $workcond/mean
# bias field correction on the mean image (this is because the t1 image is also bias field corrected and we are matching it to that)
N4BiasFieldCorrection -i $workcond/mean.nii.gz -o $workcond/bmean.nii.gz
# mask participant's T1-weighted structural image with participant's T1-space brain mask
fslmaths $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_desc-preproc_T1w.nii.gz -mas $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_desc-brain_mask.nii.gz $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-01_desc-brain_T1w.nii.gz
# coregister participant's T1 (brain only) to participant's mean functional image - rigid-body only (-t r). This generates tmp0GenericAffine.mat.
/group/mlr-lab/AH/Projects/toolboxes/ANTs/ANTs/Scripts/antsRegistrationSyN.sh -d 3 -t r -f $workcond/bmean.nii.gz -m $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-01_desc-brain_T1w.nii.gz -o $workcond/tmp -n 16 -j 1
# resample participant's T1 IN MNI SPACE (brain only) to EPI resolution (2.5mm isotropic)
flirt -in $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_*MNI*preproc_T1w.nii.gz -ref $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_*MNI*preproc_T1w.nii.gz -applyisoxfm 2.5 -interp nearestneighbour -out $workcond/tmpT1

# setup output directory
mkdir -p $dirp/derivatives/halaiprep/sub-"$ids"/func/ME_report/"$data"
# transform optimally combined (rec-t2star) functional data from native EPI space to MNI space (using resampled T1 in MNI space as a reference). First use tmp0GenericAffine.mat (inverse version), generated in the antsRegistration call above, to transform the EPI images into T1 space. Then use the .h5 files provided by fmriprep (during anatomical preprocessing) to transform into MNI space
antsApplyTransforms -d 3 -e 3 -i $workcond/udesc-optcom_bold.nii.gz -r $workcond/tmpT1.nii.gz -o $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-t2star_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz --default-value 0 --float 1 -n LanczosWindowedSinc --transform $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 --transform [$workcond/tmp0GenericAffine.mat,1] --transform identity --transform identity
# do the same with the denoised (rec-tedana data)
antsApplyTransforms -d 3 -e 3 -i $workcond/udesc-optcomDenoised_bold.nii.gz -r $workcond/tmpT1.nii.gz -o $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-tedana_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz --default-value 0 --float 1 -n LanczosWindowedSinc --transform $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5 --transform [$workcond/tmp0GenericAffine.mat,1] --transform identity --transform identity
# resample participant's brain mask IN MNI SPACE to EPI resolution (2.5mm isotropic)
flirt -in $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-1_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz -ref $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_*MNI*preproc_T1w.nii.gz -applyisoxfm 2.5 -interp nearestneighbour -out $workcond/tmpT1_mask
# mask optimally combined (rec-t2star) data
fslmaths $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-tedana_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz -mas $workcond/tmpT1_mask.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-tedana_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz
# mask denoised (rec-tedana) data
fslmaths $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-t2star_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz -mas $workcond/tmpT1_mask.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-t2star_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz

#tidy files
# copy motion text file into /derivatives
cp $workcond/*motion.txt $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_motion.txt
# copy affine transformation between native T1 space and native EPI space into /derivatives (you can use this to transform in both directions - T1 to EPI and EPI to T1)
cp $workcond/tmp0GenericAffine.mat $dirp/derivatives/fmriprep/sub-"$ids"/anat/sub-"$ids"_run-01_from-T1w_to-"$data"_mode-image_xfm.mat
# copy contents of /tedana folder into /derivatives
cp -rf $workcond/tedana/* $dirp/derivatives/halaiprep/sub-"$ids"/func/ME_report/"$data"/

done



