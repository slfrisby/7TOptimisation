function [] = mp2rage_proc(rootdir,ids)

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771/')

subcode{1} = ['sub-',ids];
s=1;

d=dir([rootdir,'/work/',subcode{s},'/',subcode{s},'*.nii']);
inv1=[rootdir,'/work/',subcode{s},'/',subcode{s},'_run-01_INV1.nii'];
inv2=[rootdir,'/work/',subcode{s},'/',subcode{s},'_run-01_INV2.nii'];
uni=[rootdir,'/work/',subcode{s},'/',subcode{s},'_run-01_UNI.nii'];

% this calls a function called RobustCombination.m which implements O'Brien regularisation (O'Brien et al. 2014) to remove background noise
removebackgroundnoise(uni,inv1,inv2)
% this simply multiplies the images together
spm_imcalc_exp(uni,inv2);

% get the filename
[apath,afile,anext]=fileparts(uni);
uni_denoised=fullfile(apath,strcat(afile,'_denoised',anext));

% make additional functionality available
cat12('expert')

matlabbatch{1}.spm.tools.cat.estwrite.data = {uni_denoised}; % specify which files to use
matlabbatch{1}.spm.tools.cat.estwrite.data_wmh = {''};
matlabbatch{1}.spm.tools.cat.estwrite.nproc = 24; % how many separate computational processes to split it into (this is slightly lower than the cat12 default of 28)
matlabbatch{1}.spm.tools.cat.estwrite.useprior = '';
matlabbatch{1}.spm.tools.cat.estwrite.opts.tpm = {'/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771/tpm/TPM.nii'}; % for initial voxel-based preprocessing, select a map that gives probabilities of grey matter, white matter etc.
matlabbatch{1}.spm.tools.cat.estwrite.opts.affreg = 'mni'; % use standard mni template for affine registration
matlabbatch{1}.spm.tools.cat.estwrite.opts.biasstr = 1; % strength of bias field correction (cat12 default 0.5; 1 is considered "heavy") 
matlabbatch{1}.spm.tools.cat.estwrite.opts.accstr = 0.5; % accuracy level demanded from spm - this is default, but you can amp it up for very inhomogenous images
% segmentation options
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.APP = 1070; % level of initial bias correction (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.NCstr = -Inf; % strength of noise correction (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.spm_kamap = 0; % default - needs adjusting only if very abnormal brains are expected
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.LASstr = 0.5; % strength of correction for inhomogeneities in grey matter intensity (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.gcutstr = 2; % method of initial skull-stripping (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.cleanupstr = 0.5; % strength of "cleanup" (removing meninges and correcting for partial volume effects) after final (AMAP) segmentation (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.BVCstr = 0.5; % blood vessel correction
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.WMHC = 3; % set white matter hyperintensity as its "own class", rather than including it in the white matter (NOT DEFAULT)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.SLC = 0; % no stroke lesion correction
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.mrf = 1; % Markov Random Field noise filtering
matlabbatch{1}.spm.tools.cat.estwrite.extopts.segmentation.restypes.optimal = [1 0.1]; % voxel size for resampling
% registration
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.shootingtpm = {'/imaging/local/software/spm_cbu_svn/releases/spm12_fil_r7771/toolbox/cat12/templates_volumes/Template_0_IXI555_MNI152_GS.nii'}; % template for shooting registration (this is just the name for the algorithm)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.registration.shooting.regstr = 0.5; % standard shooting method
% surface options
matlabbatch{1}.spm.tools.cat.estwrite.extopts.vox = 1; % voxel size for final image (NOT DEFAULT)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtres = 0.5; % default options for cortical thickness estimation and reconstruction of the central surface using a projection-based thickness method
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtmethod = 'pbt2x';
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.pbtlas = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.collcorr = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.reduce_mesh = 1;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.vdist = 1.33333333333333;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.scale_cortex = 0.7; % avoid "glued" gyri and sulci (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.add_parahipp = 0.1; % avoid large cuts in the parahippocampal gyrus (default)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.surface.close_parahipp = 0;
% admin options (all defaults)
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.experimental = 0; % just use the tried-and-tested version
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.new_release = 0;
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.lazy = 0; % redo processing even if a result already exists
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.ignoreErrors = 1; % if an error occurs, just keep going with the next dataset
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.verb = 2; % print out the details as you go along
matlabbatch{1}.spm.tools.cat.estwrite.extopts.admin.print = 0; % don't create final CAT report since you can't read it anyway
% writing options 
matlabbatch{1}.spm.tools.cat.estwrite.output.surface = 0; % don't write surface data
matlabbatch{1}.spm.tools.cat.estwrite.output.surf_measures = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ROImenu.noROI = struct([]); % don't output ROIs
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.native = 1; % save native space grey matter
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.warped = 1; % save warped (but not modulated) grey matter
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.mod = 1; % save warped and modulated grey matter
matlabbatch{1}.spm.tools.cat.estwrite.output.GM.dartel = 1; % export data for use with DARTEL rigid body transformations 
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.native = 1; % do the same for white matter
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.WM.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.native = 1; % do the same for CSF
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.warped = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.mod = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.CSF.dartel = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.native = 0; % but don't bother with percentage position, white matter hyperintensities, stroke lesions etc.
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.ct.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.pp.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.WMH.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.SL.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.native = 1; % do save a tisue probability map in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.mod = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.TPMC.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.atlas.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.label.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.labelnative = 1;
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.native = 1; % and an image corrected for bias and global inhomogeneity in native space
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.warped = 1; % and warped
matlabbatch{1}.spm.tools.cat.estwrite.output.bias.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.native = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.warped = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.las.dartel = 0;
matlabbatch{1}.spm.tools.cat.estwrite.output.jacobianwarped = 0; 
matlabbatch{1}.spm.tools.cat.estwrite.output.warps = [0 0]; 
matlabbatch{1}.spm.tools.cat.estwrite.output.rmat = 0;

spm_jobman('run',matlabbatch);

% create and apply brain mask. We don't use this masked version, but it is
% useful as a sanity check - we can compare it to the output of fMRIprep
% and check that fMRIprep really does a better job.

% get results of segmentation - first the ones that aren't modulated or
% normalised
d=dir([rootdir,'/work/',subcode{s},'/mri/p*.nii']);
% get all components except the air outside the head, plus the modulated
% and normalised components
d(end:end+3,1)=dir([rootdir,'/work/',subcode{s},'/mri/m*.nii']);
% combine to create bias-corrected brain
Vo = spm_imcalc(spm_vol([[d(1).folder,filesep,d(2).name];[d(1).folder,filesep,d(3).name];[d(1).folder,filesep,d(4).name]]),[d(1).folder,filesep,'coreg_filled_brainmask.nii'],'imfill(i1+i2+i3>0,''holes'')');
% combine to create brain mask
Vo = spm_imcalc([spm_vol([d(1).folder,filesep,'coreg_filled_brainmask.nii']);spm_vol([d(1).folder,filesep,d(7).name])],[d(1).folder,filesep,'bias_corrected_brain.nii'],'imfill(i1>0,''holes'').*i2');

end
        
