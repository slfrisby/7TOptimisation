addpath('/group/mlr-lab/AH/Projects/spm12/');
addpath(genpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/riksneurotools-master/'));

root='/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives';

% exclusions: 013 and 017 for excessive head motion, and 004 because the
% pTx run was not acquired.
subs=[{'001'},{'002'},{'003'},{'005'},{'006'},{'007'},{'008'},{'009'},{'010'},{'011'},{'012'},{'014'},{'015'},{'016'},{'017'},{'018'},{'019'}];

%% comparing pTx to SESB

% model fit (con images)

% setup files for batch_spm_anova. This setup is slightly different from
% the recommended setup in the comments for batch _spm_anova, but produces
% the correct design matrix.
imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_ptx8ms/con_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/con/SESB_vs_pTx/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'pTx';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1];
S.contrasts{n}.name = 'SESB>pTx';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1]; 
S.contrasts{n}.name = 'pTx>SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,2) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

% model precision (spmT images)

% setup files for batch_spm_anova.
imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_ptx8ms/spmT_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/spmT/SESB_vs_pTx/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'pTx';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1];
S.contrasts{n}.name = 'SESB>pTx';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1]; 
S.contrasts{n}.name = 'pTx>SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,2) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);


%% 2x2 Factorial design - echo and band

% model fit (con images)

% Setup files for batch_spm_anova.
imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB/con_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/con_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/con/factorial/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 1 -1 -1];
S.contrasts{n}.name = 'SE>ME';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 -1 1 1]; 
S.contrasts{n}.name = 'ME>SE';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'SB>MB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1]; 
S.contrasts{n}.name = 'MB>SB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

% model precision (spmT images)

imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB/spmT_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/spmT_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/spmT/factorial/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 1 -1 -1];
S.contrasts{n}.name = 'SE>ME';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 -1 1 1]; 
S.contrasts{n}.name = 'ME>SE';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'SB>MB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1]; 
S.contrasts{n}.name = 'MB>SB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

%% Effect of denoising

% model fit (con images)

% Setup files for batch_spm_anova.
imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB_dn/con_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_dn/con_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/con/denoising/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB_dn';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'MEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MEMB_dn';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'standard>denoised';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1];
S.contrasts{n}.name = 'denoised>standard';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

% model precision (spmT images)

imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB_dn/spmT_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_dn/spmT_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/spmT/denoising/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MESB_dn';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'MEMB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MEMB_dn';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'standard>denoised';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1];
S.contrasts{n}.name = 'denoised>standard';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

%% Odd volumes

% model fit (con images)

% Setup files for batch_spm_anova.
imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB_odd/con_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_odd/con_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/con/odd_volumes/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SEMB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'MEMB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 1 -1 -1]; 
S.contrasts{n}.name = 'SE>ME';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 -1 1 1]; 
S.contrasts{n}.name = 'ME>SE';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'SB>MB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1]; 
S.contrasts{n}.name = 'MB_odd>SB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);

% model precision (spmT images)

imgs = cell(1,length(subs));
S = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB_odd/spmT_0003.nii']; 
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_odd/spmT_0003.nii']; 
 S.imgfiles{1}{s} = strvcat(imgs{s});
end

S.outdir = [root,'/GLM/second/spmT/odd_volumes/S_gt_C'];
mkdir(S.outdir);
S.mask = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii'];

n=1;
S.contrasts{n}.c = [1 0 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 1 0 0 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'SEMB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 1 0 ones(1,length(subs))/length(subs)];
S.contrasts{n}.name = 'MESB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [0 0 0 1 ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'MEMB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 1 -1 -1]; 
S.contrasts{n}.name = 'SE>ME';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 -1 1 1]; 
S.contrasts{n}.name = 'ME>SE';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 1 -1]; 
S.contrasts{n}.name = 'SB>MB_odd';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [-1 1 -1 1]; 
S.contrasts{n}.name = 'MB_odd>SB';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [1 -1 -1 1]; 
S.contrasts{n}.name = 'Interaction';
S.contrasts{n}.type = 'T';
n=n+1;
S.contrasts{n}.c = [ones(1,4) ones(1,length(subs))/length(subs)]; 
S.contrasts{n}.name = 'maineffect';
S.contrasts{n}.type = 'T';
n=n+1;

spm('defaults','fMRI');
batch_spm_anova(S);
