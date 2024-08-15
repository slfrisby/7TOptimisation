#!/bin/bash
echo "
++++++++++++++++++++++++" 

#ids is set in sub_job.sh 
#ids=$1
dirp=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main
work=/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work

# runs 1st level GLMs for main analysis
#matlab_r2019a -nodisplay -nodesktop -r "addpath('$dirp/scripts/');firstlevel_glm('"$ids"',$sm);exit"
# runs 1st level GLMs for MVPA
#matlab_r2019a -nodisplay -nodesktop -r "addpath('$dirp/scripts/');firstlevel_glm_mvpa('"$ids"',$sm);exit"
# runs 1st level GLMs for slice leakage artifact testing
matlab_r2019a -nodisplay -nodesktop -r "addpath('$dirp/scripts/');firstlevel_glm_native('"$ids"');exit"


