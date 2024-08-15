function [] = slice_leakage_make_volume_mask(workdir)
% make a "volume mask" (matrix the same size as a brain volume, all filled
% with 1s) for use in t2star combination during preprocessing for slice
% leakage testing.

mask = ones(84,84,48);
niftiwrite(mask,[workdir,'/volume_mask.nii']);
gzip([workdir,'/volume_mask.nii']);

end