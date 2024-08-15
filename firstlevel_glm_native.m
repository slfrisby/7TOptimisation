function [subcode] = firstlevel_glm_native(ids)

addpath('/group/mlr-lab/AH/Projects/spm12/')
root = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/'];
cd([root]);

subcode{1} = ['sub-',ids] 

%% run 1st level analysis for all subjects

cond={'SESB','SEMB','MESB','MEMB'};
tr=[3.02,1.51,3.02,1.51];
nslices=[48,48,48,48,48];

% for every participant
for s=1:size(subcode,2)
    
    % for every sequence
    for c=1:size(cond,2)
        cd(root)

        % find motion parameter file and set directories
        if exist([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt'])
            disp([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt','*****FILE FOUND - PROCESSING*****'])
            outdir=([root,'GLM/first_native/',subcode{s},'/',cond{c},'/'])
            datadir=([root,'SPM/first_native/',subcode{s}]);
            mkdir(outdir);
            mkdir(datadir);
            delete([outdir,'SPM.mat']);

            if isempty(dir([datadir,'/',cond{c},'*','.nii']))
                %load confounds file and extract 6 motion parameters
                R=spm_load(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt']);
                save([datadir,'/motion_',cond{c},'.mat'],'R');  

                %unzip func file to be used in SPM
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'*_space-native_*_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name); % this gets the UNZIPPED name (tmp(1).name gives the zipped name)
                gunzip([tmp(1).folder,'/',tmp(1).name],[datadir,'/']);
            end

            %select design matrix
            des_mat=([root,'../scripts/design_matrix.mat']);  %need to create with original eprime output once

            %build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr(c); % TR variable set above
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nslices(c); % slices variable set above
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = nslices(c)/2; % align predictor variables so they are predicting responses midway through volume acquisition
            switch c % specify scans
                case 1
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:171])); %load unsmoothed data. Number of volumes change by protocol but are consistent between participants
                case 2
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:340]));   
                case 3
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:171])); 
                case 4
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:340]));      
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'.mat']}; % specify motion regressors - they are extracted from the confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
            % no mask
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian
            matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'S';
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'C';
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [0 1];
            matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.name = 'S>C';
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
            matlabbatch{3}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'C>S';
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
            matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'S+C';
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [1 1];
            matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.delete = 1;
            spm_jobman('run',matlabbatch);
            
            % also build GLM for MVPA slice leakage testing
            outdir=([root,'GLM/first_native_mvpa/',subcode{s},'/',cond{c},'/'])
            datadir=([root,'SPM/first_native_mvpa/',subcode{s}]);
            mkdir(outdir);
            mkdir(datadir);
            delete([outdir,'SPM.mat']);
            
            if isempty(dir([datadir,'/',cond{c},'*','.nii']))
                %load confounds file and extract 6 motion parameters
                R=spm_load(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt']);
                save([datadir,'/motion_',cond{c},'.mat'],'R');  

                %unzip func file to be used in SPM
                tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'*_space-native_*_bold.nii.gz']);
                [~,name,~]=fileparts(tmp(1).name); % this gets the UNZIPPED name (tmp(1).name gives the zipped name)
                gunzip([tmp(1).folder,'/',tmp(1).name],[datadir,'/']);
            end
            
            des_mat=([root,'../scripts/design_matrix_mvpa.mat']);  %need to create with original eprime output once

            %build GLM
            clear matlabbatch
            matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
            matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
            matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr(c); % TR variable set above
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nslices(c); % slices variable set above
            matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = nslices(c)/2; % align predictor variables so they are predicting responses midway through volume acquisition
            switch c % specify scans
                case 1
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:171])); %load unsmoothed data. Number of volumes change by protocol but are consistent between participants
                case 2
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:340]));   
                case 3
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:171])); 
                case 4
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],[name],[1:340]));      
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'.mat']}; % specify motion regressors - they are extracted from the confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
            % no mask
            matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; % high-pass filter - slow signal drifts with a period LONGER than this will be removed (this is the default value)
            matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {}); % format to expect design matrix in
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {}); % format to expect motion regressors in
            matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
            matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0]; % do not assume that the HRF varies between subjects and voxels
            matlabbatch{1}.spm.stats.fmri_spec.volt = 1; % use linear convolution for the HRF (https://pubmed.ncbi.nlm.nih.gov/9438436/)
            matlabbatch{1}.spm.stats.fmri_spec.global = 'None'; % no global intensity normalisation
            matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0; % default to go with the above
            matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)'; % model autocorrelation

            matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1; % as opposed to Bayesian
            % no contrasts
            spm_jobman('run',matlabbatch);

        end
    end
end

