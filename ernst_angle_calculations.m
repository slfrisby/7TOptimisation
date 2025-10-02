addpath('/group/mlr-lab/AH/Projects/spm12/');
addpath('/group/mlr-lab/AH/Projects/toolboxes/');
addpath('/group/mlr-lab/AH/Projects/toolboxes/Violinplot/');
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/riksneurotools-master/Util');
% the above line adds the function roi_extract.m to path. Before using this
% % function, check that lines 195-196 read:
%                 ROI(nr,s).mean   = nanmean(d,2);
%                 ROI(nr,s).median = nanmedian(d,2);
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts');

root = ['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives'];
cd(root);

% setup ROIs. Use all ROIs that overlapped with any main effect mentioned
% in the main text
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{7}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri

%% sub-021

% extract ROIs from con images (expressing model fit)

con = struct();

imgs{1}{1} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SESB/con_0003.nii']; 
imgs{1}{2} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SESBernst/con_0003.nii'];
imgs{1}{3} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SEMB/con_0003.nii']; 
imgs{1}{4} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SEMBernst/con_0003.nii']; 

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

spmT = struct();

imgs{1}{1} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SESB/spmT_0003.nii']; 
imgs{1}{2} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SESBernst/spmT_0003.nii'];
imgs{1}{3} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SEMB/spmT_0003.nii']; 
imgs{1}{4} = [root,'/GLM/first_mni_ernst/sub-021/6sm_SEMBernst/spmT_0003.nii']; 

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

% extract ROIs from tSNR images 

tSNR = struct();

imgs{1}{1} = [root,'/tSNR/sub-021/SESB_tSNR.nii.gz']; 
imgs{1}{2} = [root,'/tSNR/sub-021/SESBernst_tSNR.nii.gz']; 
imgs{1}{3} = [root,'/tSNR/sub-021/SEMB_tSNR.nii.gz']; 
imgs{1}{4} = [root,'/tSNR/sub-021/SEMBernst_tSNR.nii.gz']; 

% extract data from each ROI
tSNR.Datafiles = imgs;
tSNR.output_raw = 1;
tSNR.ROIfiles=R.ROIfiles;
tSNR.ROI = roi_extract(tSNR);


%remove NaN voxels from data
% for each protocol
for c=1:length(imgs{1,1})
    % for each ROI
    for n=1:length(con.ROIfiles)
        % find NaNs in the contrast
        x=isnan(con.ROI(n).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data and the corresponding values
        % from the spmT data
        con.ROI(n).rawdata(:,ind)=[];
        spmT.ROI(n).rawdata(:,ind)=[];
        tSNR.ROI(n).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con.ROI(n).XYZ(:,ind)=[];
        spmT.ROI(n).XYZ(:,ind)=[];
        tSNR.ROI(n).XYZ(:,ind)=[];
        % extract the mean and median for each ROI (for each contrast)
        con_collate_median{n}(1,:)=con.ROI(n).median';
        con_collate_mean{n}(1,:)=con.ROI(n).mean';
        spmT_collate_median{n}(1,:)=spmT.ROI(n).median';
        spmT_collate_mean{n}(1,:)=spmT.ROI(n).mean';
        tSNR_collate_median{n}(1,:)=tSNR.ROI(n).median';
        tSNR_collate_mean{n}(1,:)=tSNR.ROI(n).mean';
    end
end

% make a table of con results

cond = {'SESB','SESBernst', 'SESB % change', 'SEMB','SEMBernst', 'SEMB % change'};

for i=1:length(R.ROIfiles)
    % fill in results
    conResults(i,[1,2,4,5])=con_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    conResults(i,3) = ((conResults(i,2) - conResults(i,1)) / abs(conResults(i,1))) * 100;
    conResults(i,6) = ((conResults(i,5) - conResults(i,4)) / abs(conResults(i,4))) * 100;
    % conResults(i,3) = ((conResults(i,1) - conResults(i,2)) / abs(conResults(i,2))) * 100;
    % conResults(i,6) = ((conResults(i,4) - conResults(i,5)) / abs(conResults(i,5))) * 100;
end
conResults=array2table(conResults,'VariableNames',cond);

% make a table of spmT results

for i=1:length(R.ROIfiles)
    % fill in results
    spmTResults(i,[1,2,4,5])=spmT_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    spmTResults(i,3) = ((spmTResults(i,2) - spmTResults(i,1)) / abs(spmTResults(i,1))) * 100;
    spmTResults(i,6) = ((spmTResults(i,5) - spmTResults(i,4)) / abs(spmTResults(i,4))) * 100;
    % spmTResults(i,3) = ((spmTResults(i,1) - spmTResults(i,2)) / abs(spmTResults(i,2))) * 100;
    % spmTResults(i,6) = ((spmTResults(i,4) - spmTResults(i,5)) / abs(spmTResults(i,5))) * 100;
end
spmTResults=array2table(spmTResults,'VariableNames',cond);

% make a table of tSNR results

for i=1:length(R.ROIfiles)
    % fill in results
    tSNRResults(i,[1,2,4,5])=tSNR_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    tSNRResults(i,3) = ((tSNRResults(i,2) - tSNRResults(i,1)) / abs(tSNRResults(i,1))) * 100;
    tSNRResults(i,6) = ((tSNRResults(i,5) - tSNRResults(i,4)) / abs(tSNRResults(i,4))) * 100;
    % tSNRResults(i,3) = ((tSNRResults(i,1) - tSNRResults(i,2)) / abs(tSNRResults(i,2))) * 100;
    % tSNRResults(i,6) = ((tSNRResults(i,4) - tSNRResults(i,5)) / abs(tSNRResults(i,5))) * 100;
end
tSNRResults=array2table(tSNRResults,'VariableNames',cond);

% extract ROIs in preparation for MVPA

cond = {'SESB','SESBernst','SEMB','SEMBernst'};
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)

    % setup beta image structure
    imgs{1}{1} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0001.nii'];
    imgs{1}{2} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0002.nii'];
    imgs{1}{3} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0003.nii'];
    imgs{1}{4} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0004.nii'];
    imgs{1}{5} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0005.nii'];
    imgs{1}{6} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0006.nii'];
    imgs{1}{7} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0007.nii'];
    imgs{1}{8} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0008.nii'];
    imgs{1}{9} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0009.nii'];
    imgs{1}{10} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0010.nii'];
    imgs{1}{11} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0011.nii'];
    imgs{1}{12} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0012.nii'];
    imgs{1}{13} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0013.nii'];
    imgs{1}{14} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0014.nii'];
    imgs{1}{15} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0015.nii'];
    imgs{1}{16} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0016.nii'];
    imgs{1}{17} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0017.nii'];
    imgs{1}{18} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0018.nii'];
    imgs{1}{19} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0019.nii'];
    imgs{1}{20} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0020.nii'];
    imgs{1}{21} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0021.nii'];
    imgs{1}{22} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0022.nii'];
    imgs{1}{23} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0023.nii'];
    imgs{1}{24} = [root,'/GLM/first_mni_ernst_mvpa/sub-021/',cond{c},'/beta_0024.nii'];


% extract data from each ROI and store in one big struct
field = cond{c};
con_mvpa.(field).Datafiles = imgs;
con_mvpa.(field).output_raw = 1;
con_mvpa.(field).ROIfiles = R.ROIfiles;
con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));

end

% remove NaN voxels from data
% for each protocol
for c=1:length(cond)
    % for each ROI
    for n=1:length(R.ROIfiles)
        field = cond{c};
        % find NaNs in the beta images
        x=isnan(con_mvpa.(field).ROI(n).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data
        con_mvpa.(field).ROI(n).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con_mvpa.(field).ROI(n).XYZ(:,ind)=[];
    end
end

% for each protocol
for c=1:length(cond)
    % for each ROI
    for n=1:length(R.ROIfiles)
        field = cond{c};
        % extract data (this matrix is beta images x nonzero voxels)
        x=con_mvpa.(field).ROI(n).rawdata;
        % Each block (12 semantic and 12 control) is one row of x.
        % Calculate cosine distance between each pair of blocks (This
        % is stored as the upper triangle of the similarity matrix). 
        con_mvpa.(field).ROI(n).dissimilarity=triu(squareform(pdist(x,'cosine')));
        % convert zeros in the matrix to NaNs
        con_mvpa.(field).ROI(n).dissimilarity(con_mvpa.(field).ROI(n).dissimilarity==0)=nan;
        % calculate the mean dissimilarity between pairs of blocks in
        % the same condition (i.e. mean dissimilarity between pairs of
        % semantic blocks or pairs of control blocks)
        con_mvpa.(field).ROI(n).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n).dissimilarity([13:24],[13:24]),[],1)]);
        % calculate the mean dissimilarity between pairs of blocks in
        % different conditions (i.e. mean dissimilarity between one
        % semantic block and one control block)
        con_mvpa.(field).ROI(n).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n).dissimilarity([1:12],[13:24]),[],1));
        % calculate the difference in means
        con_mvpa.(field).ROI(n).mvpa_comparison_mean=[con_mvpa.(field).ROI(n).mvpa_between_mean-con_mvpa.(field).ROI(n).mvpa_within_mean];
        % collate results
        con_mvpa_collate_mean{n}(c)=con_mvpa.(field).ROI(n).mvpa_comparison_mean;
    end
end

% make a table of MVPA results

cond = {'SESB','SESBernst', 'SESB % change', 'SEMB','SEMBernst', 'SEMB % change'};

for i=1:length(R.ROIfiles)
    % fill in results
    MVPAResults(i,[1,2,4,5])=con_mvpa_collate_mean{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    MVPAResults(i,3) = ((MVPAResults(i,2) - MVPAResults(i,1)) / abs(MVPAResults(i,1))) * 100;
    MVPAResults(i,6) = ((MVPAResults(i,5) - MVPAResults(i,4)) / abs(MVPAResults(i,4))) * 100;
    % MVPAResults(i,3) = ((MVPAResults(i,1) - MVPAResults(i,2)) / abs(MVPAResults(i,2))) * 100;
    % MVPAResults(i,6) = ((MVPAResults(i,4) - MVPAResults(i,5)) / abs(MVPAResults(i,5))) * 100;
end
MVPAResults=array2table(MVPAResults,'VariableNames',cond);


clear imgs conResults spmTResults tSNRResults MVPAResults con_collate_mean con_collate_median spmT_collate_mean spmT_collate_median tSNR_collate_mean tSNR_collate_median con_mvpa_collate_mean 

%% sub-022

% extract ROIs from con images (expressing model fit)

con = struct();

imgs{1}{1} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SESB/con_0003.nii']; 
imgs{1}{2} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SESBernst/con_0003.nii'];
imgs{1}{3} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SEMB/con_0003.nii']; 
imgs{1}{4} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SEMBernst/con_0003.nii']; 

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

spmT = struct();

imgs{1}{1} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SESB/spmT_0003.nii']; 
imgs{1}{2} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SESBernst/spmT_0003.nii'];
imgs{1}{3} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SEMB/spmT_0003.nii']; 
imgs{1}{4} = [root,'/GLM/first_mni_ernst/sub-022/6sm_SEMBernst/spmT_0003.nii']; 

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

% extract ROIs from tSNR images 

tSNR = struct();

imgs{1}{1} = [root,'/tSNR/sub-022/SESB_tSNR.nii.gz']; 
imgs{1}{2} = [root,'/tSNR/sub-022/SESBernst_tSNR.nii.gz']; 
imgs{1}{3} = [root,'/tSNR/sub-022/SEMB_tSNR.nii.gz']; 
imgs{1}{4} = [root,'/tSNR/sub-022/SEMBernst_tSNR.nii.gz']; 

% extract data from each ROI
tSNR.Datafiles = imgs;
tSNR.output_raw = 1;
tSNR.ROIfiles=R.ROIfiles;
tSNR.ROI = roi_extract(tSNR);

%remove NaN voxels from data
% for each protocol
for c=1:length(imgs{1,1})
    % for each ROI
    for n=1:length(con.ROIfiles)
        % find NaNs in the contrast
        x=isnan(con.ROI(n).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data and the corresponding values
        % from the spmT data
        con.ROI(n).rawdata(:,ind)=[];
        spmT.ROI(n).rawdata(:,ind)=[];
        tSNR.ROI(n).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con.ROI(n).XYZ(:,ind)=[];
        spmT.ROI(n).XYZ(:,ind)=[];
        tSNR.ROI(n).XYZ(:,ind)=[];
        % extract the mean and median for each ROI (for each contrast)
        con_collate_median{n}(1,:)=con.ROI(n).median';
        con_collate_mean{n}(1,:)=con.ROI(n).mean';
        spmT_collate_median{n}(1,:)=spmT.ROI(n).median';
        spmT_collate_mean{n}(1,:)=spmT.ROI(n).mean';
        tSNR_collate_median{n}(1,:)=tSNR.ROI(n).median';
        tSNR_collate_mean{n}(1,:)=tSNR.ROI(n).mean';
    end
end

% make a table of con results

cond = {'SESB','SESBernst', 'SESB % change', 'SEMB','SEMBernst', 'SEMB % change'};

for i=1:length(R.ROIfiles)
    % fill in results
    conResults(i,[1,2,4,5])=con_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    conResults(i,3) = ((conResults(i,2) - conResults(i,1)) / abs(conResults(i,1))) * 100;
    conResults(i,6) = ((conResults(i,5) - conResults(i,4)) / abs(conResults(i,4))) * 100;
    % conResults(i,3) = ((conResults(i,1) - conResults(i,2)) / abs(conResults(i,2))) * 100;
    % conResults(i,6) = ((conResults(i,4) - conResults(i,5)) / abs(conResults(i,5))) * 100;
end
conResults=array2table(conResults,'VariableNames',cond);

% make a table of spmT results

for i=1:length(R.ROIfiles)
    % fill in results
    spmTResults(i,[1,2,4,5])=spmT_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    spmTResults(i,3) = ((spmTResults(i,2) - spmTResults(i,1)) / abs(spmTResults(i,1))) * 100;
    spmTResults(i,6) = ((spmTResults(i,5) - spmTResults(i,4)) / abs(spmTResults(i,4))) * 100;
    % spmTResults(i,3) = ((spmTResults(i,1) - spmTResults(i,2)) / abs(spmTResults(i,2))) * 100;
    % spmTResults(i,6) = ((spmTResults(i,4) - spmTResults(i,5)) / abs(spmTResults(i,5))) * 100;
end

for i=1:length(R.ROIfiles)
    % fill in results
    tSNRResults(i,[1,2,4,5])=tSNR_collate_median{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    tSNRResults(i,3) = ((tSNRResults(i,2) - tSNRResults(i,1)) / abs(tSNRResults(i,1))) * 100;
    tSNRResults(i,6) = ((tSNRResults(i,5) - tSNRResults(i,4)) / abs(tSNRResults(i,4))) * 100;
    % tSNRResults(i,3) = ((tSNRResults(i,1) - tSNRResults(i,2)) / abs(tSNRResults(i,2))) * 100;
    % tSNRResults(i,6) = ((tSNRResults(i,4) - tSNRResults(i,5)) / abs(tSNRResults(i,5))) * 100;
end

% extract ROIs in preparation for MVPA

cond = {'SESB','SESBernst','SEMB','SEMBernst'};
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)

    % setup beta image structure
    imgs{1}{1} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0001.nii'];
    imgs{1}{2} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0002.nii'];
    imgs{1}{3} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0003.nii'];
    imgs{1}{4} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0004.nii'];
    imgs{1}{5} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0005.nii'];
    imgs{1}{6} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0006.nii'];
    imgs{1}{7} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0007.nii'];
    imgs{1}{8} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0008.nii'];
    imgs{1}{9} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0009.nii'];
    imgs{1}{10} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0010.nii'];
    imgs{1}{11} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0011.nii'];
    imgs{1}{12} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0012.nii'];
    imgs{1}{13} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0013.nii'];
    imgs{1}{14} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0014.nii'];
    imgs{1}{15} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0015.nii'];
    imgs{1}{16} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0016.nii'];
    imgs{1}{17} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0017.nii'];
    imgs{1}{18} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0018.nii'];
    imgs{1}{19} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0019.nii'];
    imgs{1}{20} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0020.nii'];
    imgs{1}{21} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0021.nii'];
    imgs{1}{22} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0022.nii'];
    imgs{1}{23} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0023.nii'];
    imgs{1}{24} = [root,'/GLM/first_mni_ernst_mvpa/sub-022/',cond{c},'/beta_0024.nii'];


% extract data from each ROI and store in one big struct
field = cond{c};
con_mvpa.(field).Datafiles = imgs;
con_mvpa.(field).output_raw = 1;
con_mvpa.(field).ROIfiles = R.ROIfiles;
con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));

end

% remove NaN voxels from data
% for each protocol
for c=1:length(cond)
    % for each ROI
    for n=1:length(R.ROIfiles)
        field = cond{c};
        % find NaNs in the beta images
        x=isnan(con_mvpa.(field).ROI(n).rawdata(1,:));
        % get indices of these NaNs
        ind=find(x);
        % remove the NaNs from the data
        con_mvpa.(field).ROI(n).rawdata(:,ind)=[];
        % remove the coordinates of the NaNs from the matrix of
        % coordinates, so that everything lines up
        con_mvpa.(field).ROI(n).XYZ(:,ind)=[];
    end
end

% for each protocol
for c=1:length(cond)
    % for each ROI
    for n=1:length(R.ROIfiles)
        field = cond{c};
        % extract data (this matrix is beta images x nonzero voxels)
        x=con_mvpa.(field).ROI(n).rawdata;
        % Each block (12 semantic and 12 control) is one row of x.
        % Calculate cosine distance between each pair of blocks (This
        % is stored as the upper triangle of the similarity matrix). 
        con_mvpa.(field).ROI(n).dissimilarity=triu(squareform(pdist(x,'cosine')));
        % convert zeros in the matrix to NaNs
        con_mvpa.(field).ROI(n).dissimilarity(con_mvpa.(field).ROI(n).dissimilarity==0)=nan;
        % calculate the mean dissimilarity between pairs of blocks in
        % the same condition (i.e. mean dissimilarity between pairs of
        % semantic blocks or pairs of control blocks)
        con_mvpa.(field).ROI(n).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n).dissimilarity([13:24],[13:24]),[],1)]);
        % calculate the mean dissimilarity between pairs of blocks in
        % different conditions (i.e. mean dissimilarity between one
        % semantic block and one control block)
        con_mvpa.(field).ROI(n).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n).dissimilarity([1:12],[13:24]),[],1));
        % calculate the difference in means
        con_mvpa.(field).ROI(n).mvpa_comparison_mean=[con_mvpa.(field).ROI(n).mvpa_between_mean-con_mvpa.(field).ROI(n).mvpa_within_mean];
        % collate results
        con_mvpa_collate_mean{n}(c)=con_mvpa.(field).ROI(n).mvpa_comparison_mean;
    end
end

% make a table of MVPA results

cond = {'SESB','SESBernst', 'SESB % change', 'SEMB','SEMBernst', 'SEMB % change'};

for i=1:length(R.ROIfiles)
    % fill in results
    MVPAResults(i,[1,2,4,5])=con_mvpa_collate_mean{1,i};  
    % calculate percentage improvement when Ernst angle used (instead of 50
    % degrees)
    MVPAResults(i,3) = ((MVPAResults(i,2) - MVPAResults(i,1)) / abs(MVPAResults(i,1))) * 100;
    MVPAResults(i,6) = ((MVPAResults(i,5) - MVPAResults(i,4)) / abs(MVPAResults(i,4))) * 100;
    % MVPAResults(i,3) = ((MVPAResults(i,1) - MVPAResults(i,2)) / abs(MVPAResults(i,2))) * 100;
    % MVPAResults(i,6) = ((MVPAResults(i,4) - MVPAResults(i,5)) / abs(MVPAResults(i,5))) * 100;
end
MVPAResults=array2table(MVPAResults,'VariableNames',cond);


clear conResults spmTResults tSNRResults con_collate_mean con_collate_median spmT_collate_mean spmT_collate_median tSNR_collate_mean tSNR_collate_median con_mvpa_collate_mean

