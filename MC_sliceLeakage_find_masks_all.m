% Takes coordinates of peak voxels and creates a spherical ROI around these
% peaks, ready for slice leakage artifact testing.
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/')
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/NIfTI_20140122')

proj_dir = '/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/';

% all participants (exclusions are conducted later)
subs={'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020'};
cond={'SEMB','MEMB'};

% for each sequence
for k=1:size(cond,2)
    
    cond{k}
    
    % for each peak
    for j = 1:5
    % load the file containing the peak coordinates for all participants.
    % MAP4SL calls this the seed, so this terminology is maintained here.
    seed_filename = [proj_dir,'/scripts/slice_leakage/voxel',num2str(j),'.mat'];
    seed=load(seed_filename);
    
        % for every participant
        for i=1:size(subs,2)

            ['sub-',subs{i}]
            
            % setup config struct for MAP4SL.
            % seed (peak) for that participant
            cfg.seed = seed.voxel(i,:);
            % radius in mm
            cfg.sphere_radius = 4;
            % construct a sphere rather than a disc
            cfg.disk = 'no';
            % default (not used in analysis)
            cfg.corr_thresh = 0.6;
            % despite the fact that we are using the same peak location for
            % every sequence, we need to make separate masks for MEMB and
            % SEMB. This is because the seed and artefact locations are
            % saved into one file. For SEMB, which has no GRAPPA, we need
            % only 2 spheres (A and B) whereas for MEMB, which has GRAPPA
            % 3, we have 4 spheres (A, A_g, B and B_g). 
            cfg.saveprefix = ['sub-' subs{i} '_' cond{k}];
            % don't flip the seed into the opposite hemisphere
            cfg.flipLR = 'no';
            % CAIPI shift - we know this from the image comment in the DICOM headers (ImageComments: 'Unaliased MB2/PE2/LB')
            cfg.Shift_FOVdevX = 2;
            % are extra controls being used? This is necessary if we find
            % evidence for slice leakage.
            cfg.extracontrolsflag = 'no';
            % set input file. This is the UNPREPROCESSED data in native
            % space. This is useful because it has a .json file of the same
            % name and MAP4SL needs to make use of this. This file is in
            % the same space as the GLM outputs used for the final
            % analyses. (This is not obvious in mricron since the header is
            % changed in preprocessing, by 3dAllineate, and so the images
            % are offset. It is evident in fslview). 
            switch cond{k}
               case 'SEMB'
                   inputfilename = [proj_dir,'data/sub-',subs{i},'/func/sub-',subs{i},'_task-semantic_acq-',cond{k},'_run-01_bold.nii.gz'];
               case 'MEMB'
                   inputfilename = [proj_dir,'data/sub-',subs{i},'/func/sub-',subs{i},'_task-semantic_acq-',cond{k},'_run-01_echo-1_bold.nii.gz'];  
            end
            % set output folder
            outputfolder=[proj_dir,'/derivatives/slice_leakage/sub-',subs{i}];
            mkdir(outputfolder);
            % run MAP4SL. Some edits have been made to this function to
            % make it compatible with WBIC scan data, to give different
            % pixel values to each ROI, and to remove
            % functionality that is not needed.
            MC_MAP4SL_extracontrols(cfg, inputfilename, outputfolder);

        end
    end
end



