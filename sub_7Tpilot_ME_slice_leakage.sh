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
$dirp/scripts/./convert_to_BIDS.sh "$ids"
#run fmriprep anat-only as func runs coreg fails randomly
$dirp/scripts/./01_anat_proc.sh "$ids"
else
echo "MP2RAGE processed files found - skipping preprocessing"
fi


#####RUN fMRI processing#####
for data in MESB MEMB; do

run=$dirp/data/sub-"$ids"/func/sub-"$ids"_task-semantic_acq-"$data"_run-01_echo-1_bold.nii.gz

filename=$(basename $run .nii.gz)

workcond="$work"/slice-leakage-sub-"$ids"/"$data"
rm -rf $workcond/*
mkdir $workcond -p

# don't despike EPIs - this command gets rid of background noise, which we want to keep for slice leakage testing. Instead, to keep the naming conventions consistent for the rest of the code, just copy the right files in and call them the right names
cp $run $workcond/tmpdata1.nii.gz
cp $(echo "${run/-1_/-2_}") $workcond/tmpdata2.nii.gz
cp $(echo "${run/-1_/-3_}") $workcond/tmpdata3.nii.gz


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

#t2star combination 
# get the value of each echo from the .json file and convert it from seconds to miliseconds
e1=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$filename.json | jq '.EchoTime')*1000 | bc)
tmp=$(echo "${filename/-1_/-2_}")
e2=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$tmp.json | jq '.EchoTime')*1000 | bc)
tmp=$(echo "${filename/-1_/-3_}")
e3=$(echo $(cat "$dirp"/data/sub-"$ids"/func/$tmp.json | jq '.EchoTime')*1000 | bc)

# No brain extraction, since we do not need a brain mask. What we do need, to make the echo combination work, is a "volume mask", which is a mask the size of the whole volume, all set to 1s
if [ ! -f $work/volume_mask.nii.gz ]; then
matlab_r2019a -nodisplay -nodesktop -r "addpath(genpath('$dirp/scripts/'));slice_leakage_make_volume_mask('"$work"');exit" 
fi

# (no tedana because, without a brain mask, the ICA will do strange things. But we do want to optimally combine the echoes, so we use this function from tedana)
t2smap -d $workcond/ratmpdata1.nii.gz $workcond/ratmpdata2.nii.gz $workcond/ratmpdata3.nii.gz -e $e1 $e2 $e3 --fittype curvefit --n-threads 16 --out-dir $workcond/tedana --mask $work/volume_mask.nii.gz
# this just makes sure that the outputs have headers in FSL's standard space. This means that transforms will be applied correctly.
fslreorient2std $workcond/tedana/desc-optcom_bold $workcond/tedana/desc-optcom_bold


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

#...and apply the field map to T2star-combined EPIs. --inindex means that the first line of acqp.text describes the main file that we want to transform (i.e. our data is acquired with AP phase encoding).
applytopup --imain=$workcond/tedana/desc-optcom_bold.nii.gz --inindex=1 --datain=$workcond/acqp.txt --topup=$workcond/tedana/topup --method=jac --out=$workcond/udesc-optcom_bold.nii.gz


#Don't coregister or normalise - if you want to see the results on a brain you should back-project the structural image into functional space. Instead, just copy the latest files over
cp $workcond/udesc-optcom_bold.nii.gz $dirp/derivatives/halaiprep/sub-"$ids"/func/sub-"$ids"_acq-"$data"_rec-t2star_run-01_space-native_desc-preproc_bold.nii.gz

done



