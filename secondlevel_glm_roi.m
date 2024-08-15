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

% exclusions: 013 and 017 for excessive head motion, and 004 because the
% pTx run was not acquired.
subs=[{'001'},{'002'},{'003'},{'005'},{'006'},{'007'},{'008'},{'009'},{'010'},{'011'},{'012'},{'014'},{'015'},{'016'},{'017'},{'018'},{'019'}];

%% comparing pTx to SESB

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from con images (expressing model fit)

imgs = cell(1,length(subs));
con = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_ptx8ms/con_0003.nii']; 
end

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

imgs = cell(1,length(subs));
spmT = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_ptx8ms/spmT_0003.nii']; 
end

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

%remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(imgs{1,1})
        % for each ROI
        for n=1:length(con.ROIfiles)
            % find NaNs in the contrast
            x=isnan(con.ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data and the corresponding values
            % from the spmT data
            con.ROI(n,s).rawdata(:,ind)=[];
            spmT.ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con.ROI(n,s).XYZ(:,ind)=[];
            spmT.ROI(n,s).XYZ(:,ind)=[];
            % extract the mean and median for each ROI (for each contrast)
            con_collate_median{n}(s,:)=con.ROI(n,s).median';
            con_collate_mean{n}(s,:)=con.ROI(n,s).mean';
            spmT_collate_median{n}(s,:)=spmT.ROI(n,s).median';
            spmT_collate_mean{n}(s,:)=spmT.ROI(n,s).mean';
        end
    end
end

% violin plot - con images

roiname=[{'LvATL'},{'RITG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'pTx'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

% for every ROI
for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,2,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i});
    ylabel('Contrast value','FontSize',20);
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark pink, pTx = purple
    v(1,1).ViolinColor = colours(7,:);
    v(1,2).ViolinColor = colours(4,:);
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);
linkaxes([ax1,ax2,ax3,ax4],'y')

% violin plots - spmT images

figure;

% for every ROI
for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=spmT_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,2,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i});
    ylabel('spmT value','FontSize',20);
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark pink, pTx = purple
    v(1,1).ViolinColor = colours(7,:);
    v(1,2).ViolinColor= colours(4,:);
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);
linkaxes([ax1,ax2,ax3,ax4],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % SESB > pTx
    % conduct paired t-test on con values
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,2)],[con_collate_median{1,n}(:,1)],'tail','left');
    % get p-values
    con_p_val(1,n)=tmp2;
    % also get t-values
    con_t_val(1,n)=tmp4.tstat;
    % repeat for spmT values
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,2)],[spmT_collate_median{1,n}(:,1)],'tail','left');
    spmT_p_val(1,n)=tmp2;
    spmT_t_val(1,n)=tmp4.tstat;
    
    % pTx > SESB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,1)],[con_collate_median{1,n}(:,2)],'tail','left');
    con_p_val(2,n)=tmp2;
    con_t_val(2,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,1)],[spmT_collate_median{1,n}(:,2)],'tail','left');
    spmT_p_val(2,n)=tmp2;
    spmT_t_val(2,n)=tmp4.tstat;
    
end

%% 2x2 Factorial design - echo and band

clear con_collate_median con_collate_mean con_p_val con_t_val spmT_collate_median spmT_collate_mean spmT_p_val spmT_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{7}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from con images (expressing model fit)

imgs = cell(1,length(subs));
con = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB/con_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/con_0003.nii']; 
end

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

imgs = cell(1,length(subs));
spmT = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB/spmT_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/spmT_0003.nii']; 
end

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

%remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(imgs{1,1})
        % for each ROI
        for n=1:length(con.ROIfiles)
            % find NaNs in the contrast
            x=isnan(con.ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data and the corresponding values
            % from the spmT data
            con.ROI(n,s).rawdata(:,ind)=[];
            spmT.ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con.ROI(n,s).XYZ(:,ind)=[];
            spmT.ROI(n,s).XYZ(:,ind)=[];
            % extract the mean and median for each ROI (for each contrast)
            con_collate_median{n}(s,:)=con.ROI(n,s).median';
            con_collate_mean{n}(s,:)=con.ROI(n,s).mean';
            spmT_collate_median{n}(s,:)=spmT.ROI(n,s).median';
            spmT_collate_mean{n}(s,:)=spmT.ROI(n,s).mean';
        end
    end
end

% violin plot - con images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LFP'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'SEMB';'MESB';'MEMB'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Contrast value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

% violin plots - spmT images

figure;

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=spmT_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('spmT value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % ME > SE
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,2)],[con_collate_median{1,n}(:,3);con_collate_median{1,n}(:,4)],'tail','left');
    % get p-values
    con_p_val(1,n)=tmp2;
    % also get t-values
    con_t_val(1,n)=tmp4.tstat;
    % repeat for spmT values
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,2)],[spmT_collate_median{1,n}(:,3);spmT_collate_median{1,n}(:,4)],'tail','left');
    spmT_p_val(1,n)=tmp2;
    spmT_t_val(1,n)=tmp4.tstat;
      
    % MB > SB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,3)],[con_collate_median{1,n}(:,2);con_collate_median{1,n}(:,4)],'tail','left');
    con_p_val(2,n)=tmp2;
    con_t_val(2,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,3)],[spmT_collate_median{1,n}(:,2);spmT_collate_median{1,n}(:,4)],'tail','left');
    spmT_p_val(2,n)=tmp2;
    spmT_t_val(2,n)=tmp4.tstat;
    
    % MEMB > MESB;
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,3)],[con_collate_median{1,n}(:,4)],'tail','left');
    con_p_val(3,n)=tmp2;
    con_t_val(3,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,3)],[spmT_collate_median{1,n}(:,4)],'tail','left');
    spmT_p_val(3,n)=tmp2;
    spmT_t_val(3,n)=tmp4.tstat;

    
end

%% Effect of denoising

clear con_collate_median con_collate_mean con_p_val con_t_val spmT_collate_median spmT_collate_mean spmT_p_val spmT_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{7}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from con images (expressing model fit)

imgs = cell(1,length(subs));
con = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB_dn/con_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_dn/con_0003.nii']; 
end

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

imgs = cell(1,length(subs));
spmT = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB_dn/spmT_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_dn/spmT_0003.nii']; 
end

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

%remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(imgs{1,1})
        % for each ROI
        for n=1:length(con.ROIfiles)
            % find NaNs in the contrast
            x=isnan(con.ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data and the corresponding values
            % from the spmT data
            con.ROI(n,s).rawdata(:,ind)=[];
            spmT.ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con.ROI(n,s).XYZ(:,ind)=[];
            spmT.ROI(n,s).XYZ(:,ind)=[];
            % extract the mean and median for each ROI (for each contrast)
            con_collate_median{n}(s,:)=con.ROI(n,s).median';
            con_collate_mean{n}(s,:)=con.ROI(n,s).mean';
            spmT_collate_median{n}(s,:)=spmT.ROI(n,s).median';
            spmT_collate_mean{n}(s,:)=spmT.ROI(n,s).mean';
        end
    end
end

% violin plot - con images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LFP'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'MESB';'MESBdn';'MEMB';'MEMBdn'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Contrast value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in
    % the paper. MESB = orange; MESBdn = yellow; MEMB = dark blue; MEMBdn =
    % pale blue
    v(1,1).ViolinColor = a(2,:);
    v(1,2).ViolinColor = a(3,:);
    v(1,3).ViolinColor = a(1,:);
    v(1,4).ViolinColor = a(6,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

% violin plots - spmT images

figure;

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=spmT_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('spmT value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. MESB = orange; MESBdn = yellow; MEMB = dark blue; MEMBdn =
    % pale blue
    v(1,1).ViolinColor = a(2,:);
    v(1,2).ViolinColor = a(3,:);
    v(1,3).ViolinColor = a(1,:);
    v(1,4).ViolinColor = a(6,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % standard > denoised
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,2);con_collate_median{1,n}(:,4)],[con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,3)],'tail','left');
    con_p_val(1,n)=tmp2;
    con_t_val(1,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,2);spmT_collate_median{1,n}(:,4)],[spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,3)],'tail','left');
    spmT_p_val(1,n)=tmp2;
    spmT_t_val(1,n)=tmp4.tstat;
    
    % denoised > standard
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,3)],[con_collate_median{1,n}(:,2);con_collate_median{1,n}(:,4)],'tail','left');
    con_p_val(2,n)=tmp2;
    con_t_val(2,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,3)],[spmT_collate_median{1,n}(:,2);spmT_collate_median{1,n}(:,4)],'tail','left');
    spmT_p_val(2,n)=tmp2;
    spmT_t_val(2,n)=tmp4.tstat;
    
end

%% Odd volumes

clear con_collate_median con_collate_mean con_p_val con_t_val spmT_collate_median spmT_collate_mean spmT_p_val spmT_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from con images (expressing model fit)

imgs = cell(1,length(subs));
con = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/con_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB_odd/con_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/con_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_odd/con_0003.nii']; 
end

% extract data from each ROI
con.Datafiles = imgs;
con.output_raw = 1;
con.ROIfiles=R.ROIfiles;
con.ROI = roi_extract(con);

% extract ROIs from spmT images (expressing model precision)

imgs = cell(1,length(subs));
spmT = struct();

for s=1:length(subs)
 imgs{s}{1} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SESB/spmT_0003.nii']; 
 imgs{s}{2} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_SEMB_odd/spmT_0003.nii'];
 imgs{s}{3} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MESB/spmT_0003.nii']; 
 imgs{s}{4} = [root,'/GLM/first_mni/sub-',subs{s},'/6sm_MEMB_odd/spmT_0003.nii']; 
end

% extract data from each ROI
spmT.Datafiles = imgs;
spmT.output_raw = 1;
spmT.ROIfiles=R.ROIfiles;
spmT.ROI = roi_extract(spmT);

%remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(imgs{1,1})
        % for each ROI
        for n=1:length(con.ROIfiles)
            % find NaNs in the contrast
            x=isnan(con.ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data and the corresponding values
            % from the spmT data
            con.ROI(n,s).rawdata(:,ind)=[];
            spmT.ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con.ROI(n,s).XYZ(:,ind)=[];
            spmT.ROI(n,s).XYZ(:,ind)=[];
            % extract the mean and median for each ROI (for each contrast)
            con_collate_median{n}(s,:)=con.ROI(n,s).median';
            con_collate_mean{n}(s,:)=con.ROI(n,s).mean';
            spmT_collate_median{n}(s,:)=spmT.ROI(n,s).median';
            spmT_collate_mean{n}(s,:)=spmT.ROI(n,s).mean';
        end
    end
end

% violin plot - con images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'SEMBodd';'MESB';'MEMBodd'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,3,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Contrast value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
ax6 = subplot(2,3,6);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

% violin plots - spmT images

figure;

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=spmT_collate_median{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,3,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('spmT value','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
ax6 = subplot(2,3,6);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % SB > MBodd
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,2);con_collate_median{1,n}(:,4)],[con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,3)],'tail','left');
    con_p_val(1,n)=tmp2;
    con_t_val(1,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,2);spmT_collate_median{1,n}(:,4)],[spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,3)],'tail','left');
    spmT_p_val(1,n)=tmp2;
    spmT_t_val(1,n)=tmp4.tstat;
    
    % MBodd > SB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_collate_median{1,n}(:,1);con_collate_median{1,n}(:,3)],[con_collate_median{1,n}(:,2);con_collate_median{1,n}(:,4)],'tail','left');
    con_p_val(2,n)=tmp2;
    con_t_val(2,n)=tmp4.tstat;
    [tmp1,tmp2,tmp3,tmp4]=ttest([spmT_collate_median{1,n}(:,1);spmT_collate_median{1,n}(:,3)],[spmT_collate_median{1,n}(:,2);spmT_collate_median{1,n}(:,4)],'tail','left');
    spmT_p_val(2,n)=tmp2;
    spmT_t_val(2,n)=tmp4.tstat;
    
end

%% Exploratory MVPA

% Comparing pTx to SESB

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from beta images

cond={'SESB';'ptx8ms'};

imgs = cell(1,length(subs));
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)
    % for each participant
    for s = 1:length(subs)
        % setup beta image structure
        imgs{s}{1} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0001.nii'];
        imgs{s}{2} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0002.nii'];
        imgs{s}{3} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0003.nii'];
        imgs{s}{4} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0004.nii'];
        imgs{s}{5} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0005.nii'];
        imgs{s}{6} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0006.nii'];
        imgs{s}{7} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0007.nii'];
        imgs{s}{8} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0008.nii'];
        imgs{s}{9} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0009.nii'];
        imgs{s}{10} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0010.nii'];
        imgs{s}{11} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0011.nii'];
        imgs{s}{12} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0012.nii'];
        imgs{s}{13} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0013.nii'];
        imgs{s}{14} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0014.nii'];
        imgs{s}{15} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0015.nii'];
        imgs{s}{16} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0016.nii'];
        imgs{s}{17} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0017.nii'];
        imgs{s}{18} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0018.nii'];
        imgs{s}{19} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0019.nii'];
        imgs{s}{20} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0020.nii'];
        imgs{s}{21} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0021.nii'];
        imgs{s}{22} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0022.nii'];
        imgs{s}{23} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0023.nii'];
        imgs{s}{24} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0024.nii'];
    
    end
    
    % extract data from each ROI and store in one big struct
    field = cond{c};
    con_mvpa.(field).Datafiles = imgs;
    con_mvpa.(field).output_raw = 1;
    con_mvpa.(field).ROIfiles = R.ROIfiles;
    con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));
    
end

% remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % find NaNs in the beta images
            x=isnan(con_mvpa.(field).ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data
            con_mvpa.(field).ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con_mvpa.(field).ROI(n,s).XYZ(:,ind)=[];
        end
    end
end

% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % extract data (this matrix is beta images x nonzero voxels)
            x=con_mvpa.(field).ROI(n,s).rawdata;
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            con_mvpa.(field).ROI(n,s).dissimilarity=triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            con_mvpa.(field).ROI(n,s).dissimilarity(con_mvpa.(field).ROI(n,s).dissimilarity==0)=nan;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            con_mvpa.(field).ROI(n,s).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n,s).dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            con_mvpa.(field).ROI(n,s).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            con_mvpa.(field).ROI(n,s).mvpa_comparison_mean=[con_mvpa.(field).ROI(n,s).mvpa_between_mean-con_mvpa.(field).ROI(n,s).mvpa_within_mean];
            % collate results
            con_mvpa_collate_mean{n}(s,c)=con_mvpa.(field).ROI(n,s).mvpa_comparison_mean;
        end
    end
end

% violin plot - beta images

roiname=[{'LvATL'},{'RITG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'pTx'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

% for every ROI
for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_mvpa_collate_mean{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,2,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i});
    ylabel('Difference in cosine distance','FontSize',20);
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark pink, pTx = purple
    v(1,1).ViolinColor = colours(7,:);
    v(1,2).ViolinColor= colours(4,:);
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,2,1);
ax2 = subplot(2,2,2);
ax3 = subplot(2,2,3);
ax4 = subplot(2,2,4);
linkaxes([ax1,ax2,ax3,ax4],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % SESB > pTx
    % conduct paired t-test on con values
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,2)],[con_mvpa_collate_mean{1,n}(:,1)],'tail','left');
    % get p-values
    con_mvpa_p_val(1,n)=tmp2;
    % also get t-values
    con_mvpa_t_val(1,n)=tmp4.tstat;

    % pTx > SESB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,1)],[con_mvpa_collate_mean{1,n}(:,2)],'tail','left');
    con_mvpa_p_val(2,n)=tmp2;
    con_mvpa_t_val(2,n)=tmp4.tstat;

end

% 2x2 factorial design - echo and band

clear con_mvpa con_mvpa_collate_mean con_mvpa_p_val con_mvpa_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{7}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from beta images

cond={'SESB';'SEMB';'MESB';'MEMB'};

imgs = cell(1,length(subs));
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)
    % for each participant
    for s = 1:length(subs)
        % setup beta image structure
        imgs{s}{1} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0001.nii'];
        imgs{s}{2} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0002.nii'];
        imgs{s}{3} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0003.nii'];
        imgs{s}{4} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0004.nii'];
        imgs{s}{5} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0005.nii'];
        imgs{s}{6} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0006.nii'];
        imgs{s}{7} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0007.nii'];
        imgs{s}{8} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0008.nii'];
        imgs{s}{9} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0009.nii'];
        imgs{s}{10} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0010.nii'];
        imgs{s}{11} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0011.nii'];
        imgs{s}{12} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0012.nii'];
        imgs{s}{13} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0013.nii'];
        imgs{s}{14} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0014.nii'];
        imgs{s}{15} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0015.nii'];
        imgs{s}{16} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0016.nii'];
        imgs{s}{17} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0017.nii'];
        imgs{s}{18} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0018.nii'];
        imgs{s}{19} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0019.nii'];
        imgs{s}{20} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0020.nii'];
        imgs{s}{21} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0021.nii'];
        imgs{s}{22} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0022.nii'];
        imgs{s}{23} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0023.nii'];
        imgs{s}{24} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0024.nii'];
    
    end
    
    % extract data from each ROI and store in one big struct
    field = cond{c};
    con_mvpa.(field).Datafiles = imgs;
    con_mvpa.(field).output_raw = 1;
    con_mvpa.(field).ROIfiles = R.ROIfiles;
    con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));
    
end

% remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % find NaNs in the beta images
            x=isnan(con_mvpa.(field).ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data
            con_mvpa.(field).ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con_mvpa.(field).ROI(n,s).XYZ(:,ind)=[];
        end
    end
end

% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % extract data (this matrix is beta images x nonzero voxels)
            x=con_mvpa.(field).ROI(n,s).rawdata;
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            con_mvpa.(field).ROI(n,s).dissimilarity=triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            con_mvpa.(field).ROI(n,s).dissimilarity(con_mvpa.(field).ROI(n,s).dissimilarity==0)=nan;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            con_mvpa.(field).ROI(n,s).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n,s).dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            con_mvpa.(field).ROI(n,s).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            con_mvpa.(field).ROI(n,s).mvpa_comparison_mean=[con_mvpa.(field).ROI(n,s).mvpa_between_mean-con_mvpa.(field).ROI(n,s).mvpa_within_mean];
            % collate results
            con_mvpa_collate_mean{n}(s,c)=con_mvpa.(field).ROI(n,s).mvpa_comparison_mean;
        end
    end
end

% violin plot - beta images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LFP'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'SEMB';'MESB';'MEMB'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_mvpa_collate_mean{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Difference in cosine distance','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % ME > SE
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,2)],[con_mvpa_collate_mean{1,n}(:,3);con_mvpa_collate_mean{1,n}(:,4)],'tail','left');
    % get p-values
    con_mvpa_p_val(1,n)=tmp2;
    % also get t-values
    con_mvpa_t_val(1,n)=tmp4.tstat;

    % MB > SB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,3)],[con_mvpa_collate_mean{1,n}(:,2);con_mvpa_collate_mean{1,n}(:,4)],'tail','left');
    con_mvpa_p_val(2,n)=tmp2;
    con_mvpa_t_val(2,n)=tmp4.tstat;
    
    % MEMB > MESB;
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,3)],[con_mvpa_collate_mean{1,n}(:,4)],'tail','left');
    con_mvpa_p_val(3,n)=tmp2;
    con_mvpa_t_val(3,n)=tmp4.tstat;

end

% Effect of denoising

clear con_mvpa con_mvpa_collate_mean con_mvpa_p_val con_mvpa_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{7}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from beta images

cond={'MESB';'MESB_dn';'MEMB';'MEMB_dn'};

imgs = cell(1,length(subs));
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)
    % for each participant
    for s = 1:length(subs)
        % setup beta image structure
        imgs{s}{1} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0001.nii'];
        imgs{s}{2} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0002.nii'];
        imgs{s}{3} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0003.nii'];
        imgs{s}{4} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0004.nii'];
        imgs{s}{5} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0005.nii'];
        imgs{s}{6} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0006.nii'];
        imgs{s}{7} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0007.nii'];
        imgs{s}{8} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0008.nii'];
        imgs{s}{9} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0009.nii'];
        imgs{s}{10} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0010.nii'];
        imgs{s}{11} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0011.nii'];
        imgs{s}{12} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0012.nii'];
        imgs{s}{13} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0013.nii'];
        imgs{s}{14} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0014.nii'];
        imgs{s}{15} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0015.nii'];
        imgs{s}{16} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0016.nii'];
        imgs{s}{17} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0017.nii'];
        imgs{s}{18} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0018.nii'];
        imgs{s}{19} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0019.nii'];
        imgs{s}{20} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0020.nii'];
        imgs{s}{21} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0021.nii'];
        imgs{s}{22} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0022.nii'];
        imgs{s}{23} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0023.nii'];
        imgs{s}{24} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0024.nii'];
    
    end
    
    % extract data from each ROI and store in one big struct
    field = cond{c};
    con_mvpa.(field).Datafiles = imgs;
    con_mvpa.(field).output_raw = 1;
    con_mvpa.(field).ROIfiles = R.ROIfiles;
    con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));
    
end

% remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % find NaNs in the beta images
            x=isnan(con_mvpa.(field).ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data
            con_mvpa.(field).ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con_mvpa.(field).ROI(n,s).XYZ(:,ind)=[];
        end
    end
end

% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % extract data (this matrix is beta images x nonzero voxels)
            x=con_mvpa.(field).ROI(n,s).rawdata;
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            con_mvpa.(field).ROI(n,s).dissimilarity=triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            con_mvpa.(field).ROI(n,s).dissimilarity(con_mvpa.(field).ROI(n,s).dissimilarity==0)=nan;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            con_mvpa.(field).ROI(n,s).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n,s).dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            con_mvpa.(field).ROI(n,s).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            con_mvpa.(field).ROI(n,s).mvpa_comparison_mean=[con_mvpa.(field).ROI(n,s).mvpa_between_mean-con_mvpa.(field).ROI(n,s).mvpa_within_mean];
            % collate results
            con_mvpa_collate_mean{n}(s,c)=con_mvpa.(field).ROI(n,s).mvpa_comparison_mean;
        end
    end
end

% violin plots - beta images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LFP'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'MESB';'MESBdn';'MEMB';'MEMBdn'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_mvpa_collate_mean{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,4,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Difference in cosine distance','FontSize',16)
    % Set violin colours, which should be consistent across all plots in
    % the paper. MESB = orange; MESBdn = yellow; MEMB = dark blue; MEMBdn =
    % pale blue
    v(1,1).ViolinColor = a(2,:);
    v(1,2).ViolinColor = a(3,:);
    v(1,3).ViolinColor = a(1,:);
    v(1,4).ViolinColor = a(6,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,4,1);
ax2 = subplot(2,4,2);
ax3 = subplot(2,4,3);
ax4 = subplot(2,4,4);
ax5 = subplot(2,4,5);
ax6 = subplot(2,4,6);
ax7 = subplot(2,4,7);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6,ax7],'y')

for n=1:length(R.ROIfiles)
    
    % standard > denoised
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,2);con_mvpa_collate_mean{1,n}(:,4)],[con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,3)],'tail','left');
    con_mvpa_p_val(1,n)=tmp2;
    con_mvpa_t_val(1,n)=tmp4.tstat;
    
    % denoised > standard
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,3)],[con_mvpa_collate_mean{1,n}(:,2);con_mvpa_collate_mean{1,n}(:,4)],'tail','left');
    con_mvpa_p_val(2,n)=tmp2;
    con_mvpa_t_val(2,n)=tmp4.tstat;
    
end

% Odd volumes

clear con_mvpa con_mvpa_collate_mean con_mvpa_p_val con_mvpa_t_val

% setup ROIs. Plot the cluster-corrected main effect (all conditions set to
% 1) and use only ROIs that overlap with this main effect (comment out
% those that do not).
R = struct();
R.ROIfiles{1}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--51_9_-36.nii'];% temporalpole
R.ROIfiles{2}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--42_-24_-27.nii']; % vATL
R.ROIfiles{3}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-45_-33_-21.nii'];% rITG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--3_57_-15.nii'];% frontal pole
R.ROIfiles{4}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--48_-21_-9.nii'];% mMTG 
R.ROIfiles{5}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--55_-43_-7.nii'];% pMTG
R.ROIfiles{6}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8--54_27_6.nii'];% IFGptri
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-51_-27_-3.nii'];% rSTG
% R.ROIfiles{n}=['/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/HumphreysPNAS2015/PNAS_semantic_sphere_8-60_0_33.nii'];%r PCG 

% extract ROIs from beta images

cond={'SESB';'SEMB_odd';'MESB';'MEMB_odd'};

imgs = cell(1,length(subs));
con_mvpa = struct();

% for each protocol
for c = 1:length(cond)
    % for each participant
    for s = 1:length(subs)
        % setup beta image structure
        imgs{s}{1} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0001.nii'];
        imgs{s}{2} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0002.nii'];
        imgs{s}{3} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0003.nii'];
        imgs{s}{4} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0004.nii'];
        imgs{s}{5} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0005.nii'];
        imgs{s}{6} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0006.nii'];
        imgs{s}{7} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0007.nii'];
        imgs{s}{8} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0008.nii'];
        imgs{s}{9} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0009.nii'];
        imgs{s}{10} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0010.nii'];
        imgs{s}{11} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0011.nii'];
        imgs{s}{12} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0012.nii'];
        imgs{s}{13} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0013.nii'];
        imgs{s}{14} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0014.nii'];
        imgs{s}{15} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0015.nii'];
        imgs{s}{16} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0016.nii'];
        imgs{s}{17} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0017.nii'];
        imgs{s}{18} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0018.nii'];
        imgs{s}{19} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0019.nii'];
        imgs{s}{20} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0020.nii'];
        imgs{s}{21} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0021.nii'];
        imgs{s}{22} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0022.nii'];
        imgs{s}{23} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0023.nii'];
        imgs{s}{24} = [root,'/GLM/first_mni_mvpa/sub-',subs{s},'/6sm_',cond{c},'/beta_0024.nii'];
    
    end
    
    % extract data from each ROI and store in one big struct
    field = cond{c};
    con_mvpa.(field).Datafiles = imgs;
    con_mvpa.(field).output_raw = 1;
    con_mvpa.(field).ROIfiles = R.ROIfiles;
    con_mvpa.(field).ROI = roi_extract(con_mvpa.(field));
    
end

% remove NaN voxels from data
% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % find NaNs in the beta images
            x=isnan(con_mvpa.(field).ROI(n,s).rawdata(1,:));
            % get indices of these NaNs
            ind=find(x);
            % remove the NaNs from the data
            con_mvpa.(field).ROI(n,s).rawdata(:,ind)=[];
            % remove the coordinates of the NaNs from the matrix of
            % coordinates, so that everything lines up
            con_mvpa.(field).ROI(n,s).XYZ(:,ind)=[];
        end
    end
end

% for each participant
for s=1:length(subs)
    % for each protocol
    for c=1:length(cond)
        % for each ROI
        for n=1:length(R.ROIfiles)
            field = cond{c};
            % extract data (this matrix is beta images x nonzero voxels)
            x=con_mvpa.(field).ROI(n,s).rawdata;
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            con_mvpa.(field).ROI(n,s).dissimilarity=triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            con_mvpa.(field).ROI(n,s).dissimilarity(con_mvpa.(field).ROI(n,s).dissimilarity==0)=nan;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            con_mvpa.(field).ROI(n,s).mvpa_within_mean=nanmean([reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[1:12]),[],1);reshape(con_mvpa.(field).ROI(n,s).dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            con_mvpa.(field).ROI(n,s).mvpa_between_mean=nanmean(reshape(con_mvpa.(field).ROI(n,s).dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            con_mvpa.(field).ROI(n,s).mvpa_comparison_mean=[con_mvpa.(field).ROI(n,s).mvpa_between_mean-con_mvpa.(field).ROI(n,s).mvpa_within_mean];
            % collate results
            con_mvpa_collate_mean{n}(s,c)=con_mvpa.(field).ROI(n,s).mvpa_comparison_mean;
        end
    end
end

% violin plot - con images

roiname=[{'LTP'},{'LvATL'},{'RITG'},{'LmMTG'},{'LpMTG'},{'LIFGpt'}];
cond={'SESB';'SEMBodd';'MESB';'MEMBodd'};
figure;

a = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for i=1:length(R.ROIfiles)
    % get data from that ROI and make it into a table
    tmp=con_mvpa_collate_mean{1,i};
    tmp=array2table(tmp,'VariableNames',cond);
    % make a violin plot
    subplot(2,3,i);
    v=violinplot(tmp);
    % set title and y-axis label
    title(roiname{1,i},'FontSize',20);
    ylabel('Difference in cosine distance','FontSize',16)
    % Set violin colours, which should be consistent across all plots in the paper. SESB = dark red; SEMB = green; MESB = orange; MEMB = dark blue
    v(1,1).ViolinColor = a(7,:);
    v(1,2).ViolinColor = a(5,:);
    v(1,3).ViolinColor = a(2,:);
    v(1,4).ViolinColor = a(1,:);
    % set axis font size
    set(gca,'FontSize',16)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
ax6 = subplot(2,3,6);
linkaxes([ax1,ax2,ax3,ax4,ax5,ax6],'y')

% statistics - post-hoc t-tests

for n=1:length(R.ROIfiles)
    
    % SB > MBodd
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,2);con_mvpa_collate_mean{1,n}(:,4)],[con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,3)],'tail','left');
    con_mvpa_p_val(1,n)=tmp2;
    con_mvpa_t_val(1,n)=tmp4.tstat;

    % MBodd > SB
    [tmp1,tmp2,tmp3,tmp4]=ttest([con_mvpa_collate_mean{1,n}(:,1);con_mvpa_collate_mean{1,n}(:,3)],[con_mvpa_collate_mean{1,n}(:,2);con_mvpa_collate_mean{1,n}(:,4)],'tail','left');
    con_mvpa_p_val(2,n)=tmp2;
    con_mvpa_t_val(2,n)=tmp4.tstat;
    
end