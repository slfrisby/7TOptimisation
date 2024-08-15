#!/bin/bash


dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
if [ ! -d $dirp/work/logs/ ]; then
mkdir -p $dirp/work/logs/
fi

#set FWHM smoothing (mm3)
sm=6

for s in 001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020; do

echo "$s"

###runs preprocessing

#runs both ME and SE workflow
#sbatch -o $dirp/work/logs/"$s"halaiprep.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/sub_7Tpilot_ME_and_SE.sh
# runs only ME workflow
#sbatch -o $dirp/work/logs/"$s"halaiprep.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/sub_7Tpilot_ME.sh
# runs only SE workflow
#sbatch -o $dirp/work/logs/"$s"halaiprep.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/sub_7Tpilot_SE.sh
# runs only anatomical preprocessing
sbatch -o $dirp/work/logs/"$s"tmp.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/01_anat_proc.sh
# runs preprocessing needed for slice leakage artifact testing
#sbatch -o $dirp/work/logs/"$s"halaiprepSL.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/sub_7Tpilot_ME_and_SE_slice_leakage.sh

### runs GLMs

#runs 1st level GLMs
#sbatch -o $dirp/work/logs/"$s"_1stglm.out -c 16 --job-name=GLM"$s" --export=ids=${s},sm=${sm} $dirp/scripts/sub_matlabjob.sh

# runs 1st level GLMs for slice leakage artifact testing
#sbatch -o $dirp/work/logs/"$s"_1stglmSL.out -c 16 --job-name=GLM"$s" --export=ids=${s} $dirp/scripts/sub_matlabjob.sh

### runs transforms

# runs transforms for slice leakage artifact testing
#sbatch -o $dirp/work/logs/"$s"transforms.out -c 16 --job-name=7T_"$s" --export=ids=${s} $dirp/scripts/slice_leakage_apply_transforms.sh

done

#001 002 003 004 005 006 007 008 009 010 011 012 013 014 015 016 017 018 019 020
	

