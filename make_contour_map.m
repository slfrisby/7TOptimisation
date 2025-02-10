% by Saskia. Creates contour maps of ROIs and whole-brain contrasts.

% N.B. This assumes that the contrasts have been plotted and saved using
% the SPM GUI. 

addpath('/group/mlr-lab/AH/Projects/spm12/')
root = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/'];
cd([root]);

% set paths to images to be converted into contours
% ROIs are, in order: temporal pole, vATL, rITG, frontal pole, mMTG, pMTG,
% IFGptri. Followed by main effect for SESB and main effect for MEMB
imgs = {[root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'],...
    [root,'/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'],...
    [root,'/derivatives/GLM/second/con/factorial/S_gt_C/SESB_0001_FWE05.nii'],...
    [root,'/derivatives/GLM/second/con/factorial/S_gt_C/MEMB_0001_FWE05.nii']};

% set the dilation kernel for dilation later
kernel = reshape([[0 1 0; 1 1 1; 0 1 0];ones(3,3);[0 1 0; 1 1 1; 0 1 0]], 3,3,3) ;

% for each image
for i = 1:length(imgs)

    % load image
    vol = spm_vol(imgs{i});
    clusters = spm_read_vols(vol);

    % binarise. ROIs are already binarised so will not be affected by this.
    % Contrast maps will be binarised
    clusters(isnan(clusters)) = 0;
    clusters(clusters~=0) = 1;

    % create contour by dilating the image (using the kernel specified above)
    % and then subtracting the original image to make a hollow shape
    contour = spm_dilate(clusters,kernel) - clusters;
    
    % set filename
    [~,tmp] = fileparts(vol.fname);
    vol.fname = [root,'/work/',tmp,'_CONTOUR.nii'];
    % set plane info
    vol.pinfo = [0;0;0];
    % save
    spm_write_vol(vol,contour);

end
