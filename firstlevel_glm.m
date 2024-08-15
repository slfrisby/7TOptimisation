function [subcode] = firstlevel_glm(ids,sm)

addpath('/group/mlr-lab/AH/Projects/spm12/')
root = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/'];
cd([root]);

subcode{1} = ['sub-',ids] 

%% run 1st level analysis for all subjects

cond={'SESB','SEMB','ptx8ms','MESB','MEMB'};
tr=[3.02,1.51,3,3.02,1.51];
nslices=[48,48,48,48,48];

% for every participant
for s=1:size(subcode,2)
    
    % for every sequence
    for c=1:size(cond,2)
        cd(root)

        % find motion parameter file and set directories
        if exist([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt'])
            disp([root,'/halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt','*****FILE FOUND - PROCESSING*****'])
            outdir=([root,'GLM/first_mni/',subcode{s},'/',num2str(sm),'sm_',cond{c},'/'])
            datadir=([root,'SPM/first_mni/',subcode{s}]);
            mkdir(outdir);
            mkdir(datadir);
            delete([outdir,'SPM.mat']);

            %if smoothed file is missing, create files
            if isempty(dir([datadir,'/s',num2str(sm),'*',cond{c},'*','.nii']))
                %load confounds file and extract 6 motion parameters
                R=spm_load(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt']);
                % for sub-001, there were some extra volumes acquired in the pTx run.
                % We want only the first 172, so we also want only the
                % first 172 values from the motion parameter file.
                if ids == '001'
                    if c == 3
                        R = R(1:172,:);
                    end
                end
                save([datadir,'/motion_',cond{c},'.mat'],'R');  

                %unzip func file to be used in SPM
                if c <= 3
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'*_bold.nii.gz']);
                else
                    %need to add this if function as ME data has two files to choose from - this part should be the standard t2star without ICA
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_rec-t2star_*_bold.nii.gz']);
                end
                [~,name,~]=fileparts(tmp(1).name); % this gets the UNZIPPED name (tmp(1).name gives the zipped name)
                gunzip([tmp(1).folder,'/',tmp(1).name],[datadir,'/']);

                % delete any old smoothed file that may exist. Wildcard prevents a
                % problem with smoothed files building up.
                if c<=3
                    delete([datadir,'/s',num2str(sm),'*',subcode{s},'_acq-',cond{c},'_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii']);
                else
                    delete([datadir,'/s',num2str(sm),'*',subcode{s},'_acq-',cond{c},'_rec-t2star_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii']);
                end
                %smooth data and delete unsmoothed
                clear matlabbatch
                matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('ExtFPList',[datadir,'/'],name,[1:999]));
                matlabbatch{1}.spm.spatial.smooth.fwhm = [sm sm sm];% in mm. This is pre-specified by when the function is called
                matlabbatch{1}.spm.spatial.smooth.dtype = 0; % specifies what sort of output file to write (0 = same as input)
                matlabbatch{1}.spm.spatial.smooth.im = 0; % don't use implicit mask (which assumes that any voxels set to 0 or NaN in the input should stay the same in the output)
                matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(sm)];
                spm_jobman('run',matlabbatch);
                delete([datadir,'/sub-*_bold.nii']);
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
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:171])); %load smoothed data. Number of volumes change by protocol but are consistent between participants
                case 2
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:340]));  
                case 3
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:172]));  
                case 4
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:171])); 
                case 5
                 matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:340]));      
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'.mat']}; % specify motion regressors - they are extracted from the confounds file
            matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
            matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'../scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii']}; % brain mask, same as template used by fMRIprep
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

            %if Multi-echo then also process the ICA-denoised datasets (doesn't do
            %anything for SE datasets)
            if c >= 4
                cd(root) 
                outdir=([root,'/GLM/first_mni/',subcode{s},'/',num2str(sm),'sm_',cond{c},'_dn/'])
                mkdir(outdir);
                delete([outdir,'/SPM.mat']);

                if isempty(dir([datadir,'/s',num2str(sm),'*',cond{c},'*tedana*.nii']))
                %unzip func file to be used in SPM
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_rec-tedana_*_bold.nii.gz']);
                    [~,name,~]=fileparts(tmp(1).name);
                    gunzip([tmp(1).folder,'/',tmp(1).name],[datadir,'/']);

                    % delete any old smoothed file that may exist
                    delete([datadir '/s',num2str(sm),'*',subcode{s},'_acq-',cond{c},'_rec-tedana_run-01_space-MNI152NLin2009cAsym_desc-preproc_bold.nii']);
                    clear matlabbatch
                    matlabbatch{1}.spm.spatial.smooth.data = cellstr(spm_select('ExtFPList',[datadir,'/'],name,[1:999]));
                    matlabbatch{1}.spm.spatial.smooth.fwhm = [sm sm sm];% in mm. This is pre-specified by when the function is called
                    matlabbatch{1}.spm.spatial.smooth.dtype = 0; % specifies what sort of output file to write (0 = same as input)
                    matlabbatch{1}.spm.spatial.smooth.im = 0; % don't use implicit mask (which assumes that any voxels set to 0 or NaN in the input should stay the same in the output)
                    matlabbatch{1}.spm.spatial.smooth.prefix = ['s',num2str(sm)];
                    spm_jobman('run',matlabbatch);
                    delete([datadir,'/sub-*_bold.nii']);
                end

                clear matlabbatch
                matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr(c); % TR variable set above
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nslices(c); % slices variable set above
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = nslices(c)/2; % align predictor variables so they are predicting responses midway through volume acquisition
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:340])); %load smoothed data 
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'.mat']}; % NOTE: motion parameters could be removed as ICA denoising supposedly removes motion related noise.
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
                matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'../scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii']}; % brain mask, same as template used by fMRIprep
                matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; % NOTE: high pass is effectively off as ICA denoising suppposedly removes frequency anomalies
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
                matlabbatch{3}.spm.stats.con.delete = 0;
                spm_jobman('run',matlabbatch);

            end

            % if multiband then also processes "odd volumes only" (for comparison
            % with single-band sequences)
            if c == 2 || c == 5
                cd(root)
                outdir=([root,'GLM/first_mni/',subcode{s},'/',num2str(sm),'sm_',cond{c},'_odd/'])
                datadir=([root,'SPM/first_mni/',subcode{s}]);
                mkdir(outdir);
                mkdir(datadir);
                delete([outdir,'SPM.mat']);
                
                % get filename
                if c == 2 
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'*_bold.nii.gz']);
                elseif c == 5
                    %need to add this if function as ME data has two files to choose from - this part should be the standard t2star without ICA
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_rec-t2star_*_bold.nii.gz']);
                end
                [~,name,~]=fileparts(tmp(1).name); % this gets the UNZIPPED name (tmp(1).name gives the zipped name)

                % extract just those motion regressors corresponding to the
                % odd volumes

                R=spm_load(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_motion.txt']);
                R = R(1:2:end,:);
                save([datadir,'/motion_',cond{c},'_odd_volumes.mat'],'R');  
    
                %select design matrix
                des_mat=([root,'../scripts/design_matrix.mat']);  %need to create with original eprime output once

                %build GLM
                clear matlabbatch
                matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
                matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr(c-1); % NOTE: modelling TR as the TR of the corresponding single-band sequence!
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nslices(c); % slices variable set above
                matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = nslices(c)/2; % align predictor variables so they are predicting responses midway through volume acquisition
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:2:340])); %load ODD VOLUMES of smoothed data
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'_odd_volumes.mat']}; % specify motion regressors - they are extracted from the confounds file
                matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
                matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'../scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii']}; % brain mask, same as template used by fMRIprep
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

                %if Multi-echo then also process the ICA-denoised datasets (doesn't do
                %anything for SE datasets)
                if c == 5
                    cd(root) 
                    outdir=([root,'/GLM/first_mni/',subcode{s},'/',num2str(sm),'sm_',cond{c},'_dn_odd/'])
                    mkdir(outdir);
                    delete([outdir,'/SPM.mat']);

                    % get filename
                    tmp=dir(['./halaiprep/',subcode{s},'/func/',subcode{s},'_acq-',cond{c},'_rec-tedana_*_bold.nii.gz']);
                    [~,name,~]=fileparts(tmp(1).name);


                    clear matlabbatch
                    matlabbatch{1}.spm.stats.fmri_spec.dir = {outdir}; % directory for outputting full design matrix (this means the final one for use in the GLM, i.e. the bit specified manually plus motion etc.)
                    matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs'; % units in which to provide timing for the parts of the design matrix that are specified manually
                    matlabbatch{1}.spm.stats.fmri_spec.timing.RT = tr(c-1); % NOTE: modelling TR as the TR of the corresponding single-band sequence!
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = nslices(c); % slices variable set above
                    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = nslices(c)/2; % align predictor variables so they are predicting responses midway through volume acquisition
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans = cellstr(spm_select('ExtFPList',[datadir,'/'],['s',num2str(sm),name],[1:2:340])); %load smoothed data 
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {[datadir,'/motion_',cond{c},'_odd_volumes.mat']}; % NOTE: motion parameters could be removed as ICA denoising supposedly removes motion related noise.
                    matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {des_mat}; % design matrix set above
                    matlabbatch{1}.spm.stats.fmri_spec.mask = {[root,'../scripts/MNI152NLin2009cAsym/tpl-MNI152NLin2009cAsym_space-MNI_res-01_brainmask.nii']}; % brain mask, same as template used by fMRIprep
                    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128; % NOTE: high pass is effectively off as ICA denoising suppposedly removes frequency anomalies
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
                    matlabbatch{3}.spm.stats.con.delete = 0;
                    spm_jobman('run',matlabbatch);


                end

            end

        end
    end
end

