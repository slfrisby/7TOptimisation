% use spherical ROIs created with MC-sliceLeakage_find_masks_all and
% extract mean contrast values within these ROIs.
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/')
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/NIfTI_20140122')
addpath('/group/mlr-lab/AH/Projects/spm12/');

proj_dir = '/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/';

% all participants (exclusions are conducted later)
subs={'001','002','003','004','005','006','007','008','009','010','011','012','013','014','015','016','017','018','019','020'};
cond={'SEMB','MEMB'};

%% mean contrast 

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
            
            % load native-space contrast images for sequence of interest
            % and control (i.e. the corresponding single-band sequence)
            seq_filename = [proj_dir,'/derivatives/GLM/first_native/sub-',subs{i},'/',cond{k},'/con_0003.nii'];
            con_filename = [proj_dir,'derivatives/GLM/first_native/sub-',subs{i},'/',cond{k}(1:2),'SB/con_0003.nii'];

            nii = load_untouch_nii(seq_filename);
            D_seq=nii.img; % get image
            nii = load_untouch_nii(con_filename);
            D_con=nii.img; % get image

            %load seed+artefact ROIs
            artefact_filename = [proj_dir,'derivatives/slice_leakage/sub-',subs{i},'/sub-',subs{i},'_',cond{k},'_seed_',num2str(seed.voxel(i,1)),'__',num2str(seed.voxel(i,2)),'__',num2str(seed.voxel(i,3)),'_sphere4_artefact_mask_conv.nii'];
            nii= load_untouch_nii(artefact_filename);
            ROI_art=nii.img;
            
            % find number of spheres. Each sphere is numbered (region A
            % contains only 1s, the first artefact region in the _label_position_link.txt file
            % contains only 2s, and so on.
            num_art = max(ROI_art,[],'all');

            for a=1:num_art
                % using nanmean since some of the contrast images contain
                % nans.
                mean_seq(i,a)=nanmean(D_seq(ROI_art==a));
                mean_con(i,a)=nanmean(D_con(ROI_art==a));
            end
        end
        
        save_filename = [proj_dir,'/derivatives/slice_leakage/',cond{k},'_seq_mean_voxel',num2str(j),'.mat'];
        save(save_filename,'mean_seq');

        save_filename = [proj_dir,'/derivatives/slice_leakage/',cond{k}(1:2),'SB_con_mean_voxel',num2str(j),'.mat'];
        save(save_filename,'mean_con');

        clear mean_seq mean_con 
    end
end

%% exploratory MVPA

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
            
            %load seed+artefact ROIs
            artefact_filename = [proj_dir,'derivatives/slice_leakage/sub-',subs{i},'/sub-',subs{i},'_',cond{k},'_seed_',num2str(seed.voxel(i,1)),'__',num2str(seed.voxel(i,2)),'__',num2str(seed.voxel(i,3)),'_sphere4_artefact_mask_conv.nii'];
            nii= load_untouch_nii(artefact_filename);
            ROI_art=nii.img;
            
            % find number of spheres. Each sphere is numbered (region A
            % contains only 1s, the first artefact region in the _label_position_link.txt file
            % contains only 2s, and so on.
            num_art = max(ROI_art,[],'all');

            % for every block
            for h = 1:24
                % load native-space contrast images for sequence of interest
                seq_filename = [proj_dir,'/derivatives/GLM/first_native_mvpa/sub-',subs{i},'/',cond{k},'/',sprintf('beta_%04d.nii',h)];
                nii = load_untouch_nii(seq_filename);
                D_seq=nii.img; 
                
                % extract beta values for every block
                for a=1:num_art
                    mvpa_data_seq(a,i).data(h,:) = D_seq(ROI_art==a)';
                end
            end
            for h = 1:24
                % load native-space contrast images for control sequence
                con_filename = [proj_dir,'/derivatives/GLM/first_native_mvpa/sub-',subs{i},'/',cond{k}(1:2),'SB/',sprintf('beta_%04d.nii',h)];
                nii = load_untouch_nii(con_filename);
                D_con=nii.img; 
                
                % extract beta values for every block
                for a=1:num_art
                    mvpa_data_con(a,i).data(h,:) = D_con(ROI_art==a)';
                end
            end

        end
        
        save_filename = [proj_dir,'/derivatives/slice_leakage/',cond{k},'_seq_mvpa_voxel',num2str(j),'.mat'];
        save(save_filename,'mvpa_data_seq');
        
        save_filename = [proj_dir,'/derivatives/slice_leakage/',cond{k}(1:2),'SB_con_mvpa_voxel',num2str(j),'.mat'];
        save(save_filename,'mvpa_data_con');

        clear mvpa_data_seq mvpa_data_con
    end
end