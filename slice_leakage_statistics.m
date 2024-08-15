% constructs violin plots and conducts t-tests to assess whether there is
% evidence for slice leakage.

addpath('/group/mlr-lab/AH/Projects/toolboxes/');
addpath('/group/mlr-lab/AH/Projects/toolboxes/Violinplot/');
% locate files
cd('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/derivatives/slice_leakage')

% EXCLUSIONS: 013 and 017 for excessive head motion, and 004 because the
% pTx run was not acquired.
subs=[{'001'},{'002'},{'003'},{'005'},{'006'},{'007'},{'008'},{'009'},{'010'},{'011'},{'012'},{'014'},{'015'},{'016'},{'017'},{'018'},{'019'}];
subsidx=logical([1,1,1,0,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,0]);

%% SEMB

% violin plots

voxelname = {'Peak 1','Peak 2', 'Peak 3', 'Peak 4', 'Peak 5'};
roiname = {'SEMB_A','SESB_A','SEMB_B','SESB_B'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    load(['SEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['SESB_con_mean_voxel',num2str(j),'.mat'])
    % reformat results into a table for plotting
    tmp = zeros(17,4);
    tmp(:,logical([1 0 1 0])) = mean_seq(subsidx,:);
    tmp(:,logical([0 1 0 1])) = mean_con(subsidx,:);
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    subplot(2,3,j);
    v = violinplot(tmp);
    % set title and y-axis label
    title(voxelname{1,j});
    ylabel('Contrast value');
    v(1,1).ViolinColor = colours(1,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(1,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(5,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(5,:);
    v(1,4).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    xtickangle(gca,-45)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
linkaxes([ax1,ax2,ax3,ax4,ax5],'y')

% statistics - t-tests

for j = 1:5
    % load sequence and artifact data
    load(['SEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['SESB_con_mean_voxel',num2str(j),'.mat'])

    % conduct t-tests
    % A
    % conduct paired t-test
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,1),mean_seq(subsidx,1),'tail','left');
    % get p-values
    SEMB_p_val(j,1) = tmp2;
    % also get t-values
    SEMB_t_val(j,1) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,2),mean_seq(subsidx,2),'tail','left');
    % get p-values
    SEMB_p_val(j,2) = tmp2;
    % also get t-values
    SEMB_t_val(j,2) = tmp4.tstat;

end

%% MEMB

% violin plots

voxelname = {'Peak 1','Peak 2', 'Peak 3', 'Peak 4', 'Peak 5'};
roiname = {'MEMB_A','MESB_A','MEMB_A_g','MESB_A_g','MEMB_B','MESB_B','MEMB_B_g','MESB_B_g'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    load(['MEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['MESB_con_mean_voxel',num2str(j),'.mat'])
    % reformat results into a table for plotting
    tmp = zeros(17,8);
    tmp(:,logical([1 0 1 0 1 0 1 0])) = mean_seq(subsidx,:);
    tmp(:,logical([0 1 0 1 0 1 0 1])) = mean_con(subsidx,:);
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    subplot(2,3,j);
    v = violinplot(tmp);
    % set title and y-axis label
    title(voxelname{1,j});
    ylabel('Contrast value');
    v(1,1).ViolinColor = colours(1,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(1,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(7,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(7,:);
    v(1,4).ViolinAlpha = 0.15;
    v(1,5).ViolinColor = colours(5,:);
    v(1,5).ViolinAlpha = 0.4;
    v(1,6).ViolinColor = colours(5,:);
    v(1,6).ViolinAlpha = 0.15;
    v(1,7).ViolinColor = colours(3,:);
    v(1,7).ViolinAlpha = 0.4;
    v(1,8).ViolinColor = colours(3,:);
    v(1,8).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    xtickangle(gca,-45)
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
linkaxes([ax1,ax2,ax3,ax4,ax5],'y')

% statistics - t-tests

for j = 1:5
    % load sequence and artifact data
    load(['MEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['MESB_con_mean_voxel',num2str(j),'.mat'])

    % conduct t-tests
    % A
    % conduct paired t-test
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,1),mean_seq(subsidx,1),'tail','left');
    % get p-values
    MEMB_p_val(j,1) = tmp2;
    % also get t-values
    MEMB_t_val(j,1) = tmp4.tstat;
    
    % A_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,2),mean_seq(subsidx,2),'tail','left');
    % get p-values
    MEMB_p_val(j,2) = tmp2;
    % also get t-values
    MEMB_t_val(j,2) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,3),mean_seq(subsidx,3),'tail','left');
    % get p-values
    MEMB_p_val(j,3) = tmp2;
    % also get t-values
    MEMB_t_val(j,3) = tmp4.tstat;
    
    % B_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,4),mean_seq(subsidx,4),'tail','left');
    % get p-values
    MEMB_p_val(j,4) = tmp2;
    % also get t-values
    MEMB_t_val(j,4) = tmp4.tstat;

end

%% exploratory MVPA

% SEMB

% set up figure
voxelname = {'Peak 1','Peak 2', 'Peak 3', 'Peak 4', 'Peak 5'};
roiname = {'SEMB_A','SESB_A','SEMB_B','SESB_B'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    
    % load data
    load(['SEMB_seq_mvpa_voxel',num2str(j),'.mat']);
    load(['SESB_con_mvpa_voxel',num2str(j),'.mat']);
    % remove excluded participants
    mvpa_data_seq = mvpa_data_seq(:,subsidx);
    mvpa_data_con = mvpa_data_con(:,subsidx);
    
    % for each participant
    for i = 1:length(subs)
        % for each artefact location
        for a = 1:2     
            
            % extract data for sequence of interest
            x = mvpa_data_seq(a,i).data;
            % remove nans
            ind = find(isnan(x(1,:)));
            x(:,ind)=[];
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            dissimilarity = triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            dissimilarity(dissimilarity==0) = NaN;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            mvpa_within_mean = nanmean([reshape(dissimilarity([1:12],[1:12]),[],1);reshape(dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            mvpa_between_mean = nanmean(reshape(dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            mvpa_comparison_mean = mvpa_between_mean - mvpa_within_mean;
            % collate results
            mvpa_collate_mean_seq(i,a) = mvpa_comparison_mean;

            % extract data for control sequence
            x = mvpa_data_con(a,i).data;
            % remove nans
            ind = find(isnan(x(1,:)));
            x(:,ind)=[];
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            dissimilarity = triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            dissimilarity(dissimilarity==0) = NaN;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            mvpa_within_mean = nanmean([reshape(dissimilarity([1:12],[1:12]),[],1);reshape(dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            mvpa_between_mean = nanmean(reshape(dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            mvpa_comparison_mean = mvpa_between_mean - mvpa_within_mean;
            % collate results
            mvpa_collate_mean_con(i,a) = mvpa_comparison_mean;
        end
    end
    
    % update violin plots
    
    % reformat results into a table for plotting
    tmp = zeros(17,4);
    tmp(:,logical([1 0 1 0])) = mvpa_collate_mean_seq;
    tmp(:,logical([0 1 0 1])) = mvpa_collate_mean_con;
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    subplot(2,3,j);
    v = violinplot(tmp);
    % set title and y-axis label
    title(voxelname{1,j});
    ylabel('Difference in cosine distance');
    v(1,1).ViolinColor = colours(1,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(1,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(5,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(5,:);
    v(1,4).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    xtickangle(gca,-45)
    
    % update t-tests
    % A
    % conduct paired t-test
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,1),mvpa_collate_mean_seq(:,1),'tail','left');
    % get p-values
    SEMB_mvpa_p_val(j,1) = tmp2;
    % also get t-values
    SEMB_mvpa_t_val(j,1) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,2),mvpa_collate_mean_seq(:,2),'tail','left');
    % get p-values
    SEMB_mvpa_p_val(j,2) = tmp2;
    % also get t-values
    SEMB_mvpa_t_val(j,2) = tmp4.tstat;
     
end

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
linkaxes([ax1,ax2,ax3,ax4,ax5],'y')

% MEMB

% set up figure
voxelname = {'Peak 1','Peak 2', 'Peak 3', 'Peak 4', 'Peak 5'};
roiname = {'MEMB_A','MESB_A','MEMB_A_g','MESB_A_g','MEMB_B','MESB_B','MEMB_B_g','MESB_B_g'};
figure;

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    
    % load data
    load(['MEMB_seq_mvpa_voxel',num2str(j),'.mat']);
    load(['MESB_con_mvpa_voxel',num2str(j),'.mat']);
    % remove excluded participants
    mvpa_data_seq = mvpa_data_seq(:,subsidx);
    mvpa_data_con = mvpa_data_con(:,subsidx);
    
    % for each participant
    for i = 1:length(subs)
        % for each artefact location
        for a = 1:4    
            
            % extract data for sequence of interest
            x = mvpa_data_seq(a,i).data;
            % remove nans
            ind = find(isnan(x(1,:)));
            x(:,ind)=[];
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            dissimilarity = triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            dissimilarity(dissimilarity==0) = NaN;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            mvpa_within_mean = nanmean([reshape(dissimilarity([1:12],[1:12]),[],1);reshape(dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            mvpa_between_mean = nanmean(reshape(dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            mvpa_comparison_mean = mvpa_between_mean - mvpa_within_mean;
            % collate results
            mvpa_collate_mean_seq(i,a) = mvpa_comparison_mean;

            % extract data for control sequence
            x = mvpa_data_con(a,i).data;
            % remove nans
            ind = find(isnan(x(1,:)));
            x(:,ind)=[];
            % Each block (12 semantic and 12 control) is one row of x.
            % Calculate cosine distance between each pair of blocks (This
            % is stored as the upper triangle of the similarity matrix). 
            dissimilarity = triu(squareform(pdist(x,'cosine')));
            % convert zeros in the matrix to NaNs
            dissimilarity(dissimilarity==0) = NaN;
            % calculate the mean dissimilarity between pairs of blocks in
            % the same condition (i.e. mean dissimilarity between pairs of
            % semantic blocks or pairs of control blocks)
            mvpa_within_mean = nanmean([reshape(dissimilarity([1:12],[1:12]),[],1);reshape(dissimilarity([13:24],[13:24]),[],1)]);
            % calculate the mean dissimilarity between pairs of blocks in
            % different conditions (i.e. mean dissimilarity between one
            % semantic block and one control block)
            mvpa_between_mean = nanmean(reshape(dissimilarity([1:12],[13:24]),[],1));
            % calculate the difference in means
            mvpa_comparison_mean = mvpa_between_mean - mvpa_within_mean;
            % collate results
            mvpa_collate_mean_con(i,a) = mvpa_comparison_mean;
        end
    end
    
    % update violin plots
    
    % reformat results into a table for plotting
    tmp = zeros(17,8);
    tmp(:,logical([1 0 1 0 1 0 1 0])) = mvpa_collate_mean_seq;
    tmp(:,logical([0 1 0 1 0 1 0 1])) = mvpa_collate_mean_con;
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    subplot(2,3,j);
    v = violinplot(tmp);
    % set title and y-axis label
    title(voxelname{1,j});
    ylabel('Difference in cosine distance');
    v(1,1).ViolinColor = colours(1,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(1,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(7,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(7,:);
    v(1,4).ViolinAlpha = 0.15;
    v(1,5).ViolinColor = colours(5,:);
    v(1,5).ViolinAlpha = 0.4;
    v(1,6).ViolinColor = colours(5,:);
    v(1,6).ViolinAlpha = 0.15;
    v(1,7).ViolinColor = colours(3,:);
    v(1,7).ViolinAlpha = 0.4;
    v(1,8).ViolinColor = colours(3,:);
    v(1,8).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    xtickangle(gca,-45)
    
    % update t-tests
    % A
    % conduct paired t-test
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,1),mvpa_collate_mean_seq(:,1),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,1) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,1) = tmp4.tstat;
    
    % A_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,2),mvpa_collate_mean_seq(:,2),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,2) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,2) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,3),mvpa_collate_mean_seq(:,3),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,3) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,3) = tmp4.tstat;
    
    % B_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,4),mvpa_collate_mean_seq(:,4),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,4) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,4) = tmp4.tstat;

end
     

% make the plots a standard size
set(gcf,'PaperUnits','centimeters','PaperSize',[20 40])
% make the axes of all subplots match
ax1 = subplot(2,3,1);
ax2 = subplot(2,3,2);
ax3 = subplot(2,3,3);
ax4 = subplot(2,3,4);
ax5 = subplot(2,3,5);
linkaxes([ax1,ax2,ax3,ax4,ax5],'y')

