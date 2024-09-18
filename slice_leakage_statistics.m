% constructs violin plots and conducts t-tests to assess whether there is
% evidence for slice leakage.
clear;
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

roiname = {'SEMB_B','SESB_B'};

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    figure;
    load(['SEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['SESB_con_mean_voxel',num2str(j),'.mat'])
    % reformat results into a table for plotting
    tmp = zeros(17,2);
    tmp(:,logical([1 0])) = mean_seq(subsidx,2);
    tmp(:,logical([0 1])) = mean_con(subsidx,2);
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    v = violinplot(tmp);
    % set y-axis label
    ylabel('Contrast value');
    v(1,1).ViolinColor = colours(5,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(5,:);
    v(1,2).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    ylim([-6,6])
    yticks(-6:2:6)
    xtickangle(gca,-45)
    set(gcf,'position',[10,10,600,600])
    saveas(gcf,['/group/mlr-lab/Saskia/7T/slice_leakage_figures/SEMB_peak',num2str(j),'_con.png'])
    close(gcf)
end

% statistics - t-tests

for j = 1:5
    % load sequence and artifact data
    load(['SEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['SESB_con_mean_voxel',num2str(j),'.mat'])
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,2),mean_seq(subsidx,2),'tail','left');
    % get p-values
    SEMB_p_val(j,1) = tmp2;
    % also get t-values
    SEMB_t_val(j,1) = tmp4.tstat;

end

%% MEMB

% violin plots

roiname = {'MEMB_A_g','SEMB_A_g','MEMB_B','MESB_B','MEMB_B_g','SEMB_B_g'};

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    figure;
    load(['MEMB_seq_mean_voxel',num2str(j),'.mat'])
    % reformat results into a table for plotting
    tmp = zeros(17,6);
    tmp(:,logical([1 0 1 0 1 0])) = mean_seq(subsidx,2:4);
    % load MESB, which is the comparison sequence for artefact region B
    load(['MESB_con_mean_voxel',num2str(j),'.mat'])
    tmp(:,logical([ 0 0 0 1 0 0])) = mean_con(subsidx,3);
    % load SEMB, which is the comparison sequence for artefact regions Ag
    % and Bg
    load(['SEMB_con_mean_voxel',num2str(j),'.mat'])
    tmp(:,logical([0 1 0 0 0 1])) = mean_congrap(subsidx,[2,4]);
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    v = violinplot(tmp);
    % set y-axis label
    ylabel('Contrast value');
    v(1,1).ViolinColor = colours(7,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(7,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(5,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(5,:);
    v(1,4).ViolinAlpha = 0.15;
    v(1,5).ViolinColor = colours(3,:);
    v(1,5).ViolinAlpha = 0.4;
    v(1,6).ViolinColor = colours(3,:);
    v(1,6).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    ylim([-6,6])
    yticks(-6:2:6)
    xtickangle(gca,-45)
    set(gcf,'position',[10,10,600,600])
    saveas(gcf,['/group/mlr-lab/Saskia/7T/slice_leakage_figures/MEMB_peak',num2str(j),'_con.png'])
    close(gcf)
end

% statistics - t-tests

for j = 1:5
    % load sequence and artifact data
    load(['MEMB_seq_mean_voxel',num2str(j),'.mat'])
    load(['MESB_con_mean_voxel',num2str(j),'.mat'])
    load(['SEMB_con_mean_voxel',num2str(j),'.mat'])

    % conduct t-tests
    
    % A_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_congrap(subsidx,2),mean_seq(subsidx,2),'tail','left');
    % get p-values
    MEMB_p_val(j,1) = tmp2;
    % also get t-values
    MEMB_t_val(j,1) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_con(subsidx,3),mean_seq(subsidx,3),'tail','left');
    % get p-values
    MEMB_p_val(j,2) = tmp2;
    % also get t-values
    MEMB_t_val(j,2) = tmp4.tstat;
    
    % B_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mean_congrap(subsidx,4),mean_seq(subsidx,4),'tail','left');
    % get p-values
    MEMB_p_val(j,3) = tmp2;
    % also get t-values
    MEMB_t_val(j,3) = tmp4.tstat;

end

%% exploratory MVPA

% SEMB

% set up figure
roiname = {'SEMB_B','SESB_B'};

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5

    figure;
    
    % load data
    load(['SEMB_seq_mvpa_voxel',num2str(j),'.mat']);
    load(['SESB_con_mvpa_voxel',num2str(j),'.mat']);
    % remove excluded participants
    mvpa_data_seq = mvpa_data_seq(:,subsidx);
    mvpa_data_con = mvpa_data_con(:,subsidx);
    
    % for each participant
    for i = 1:length(subs)
      
        % extract data for sequence of interest
        x = mvpa_data_seq(2,i).data;
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
        mvpa_collate_mean_seq(i) = mvpa_comparison_mean;

        % extract data for control sequence
        x = mvpa_data_con(2,i).data;
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
        mvpa_collate_mean_con(i) = mvpa_comparison_mean;

    end
    
    % update violin plots
    
    % reformat results into a table for plotting
    tmp = zeros(17,2);
    tmp(:,logical([1 0])) = mvpa_collate_mean_seq;
    tmp(:,logical([0 1])) = mvpa_collate_mean_con;
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    v = violinplot(tmp);
    % set title and y-axis label
    ylabel('Difference in cosine distance');
    v(1,1).ViolinColor = colours(5,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(5,:);
    v(1,2).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    ylim([-0.1,0.5])
    yticks(-0.1:0.1:0.5)
    xtickangle(gca,-45)
    set(gcf,'position',[10,10,600,600])
    saveas(gcf,['/group/mlr-lab/Saskia/7T/slice_leakage_figures/SEMB_peak',num2str(j),'_mvpa.png'])
    close(gcf)
    
    % update t-tests
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con,mvpa_collate_mean_seq,'tail','left');
    % get p-values
    SEMB_mvpa_p_val(j,1) = tmp2;
    % also get t-values
    SEMB_mvpa_t_val(j,1) = tmp4.tstat;
     
end

clear mvpa_collate_mean_seq mvpa_collate_mean_con

% MEMB

% set up figure
roiname = {'MEMB_A_g','SEMB_A_g','MEMB_B','MESB_B','MEMB_B_g','SEMB_B_g'};

colours = get(gca,'colororder'); % colours: 1 = dark blue, 2 = orange, 3 = yellow, 4 = purple, 5 = green, 6 = pale blue, 7 = dark pink/red

for j = 1:5
    figure; 

    % load data
    load(['MEMB_seq_mvpa_voxel',num2str(j),'.mat']);
    load(['MESB_con_mvpa_voxel',num2str(j),'.mat']);
    load(['SEMB_con_mvpa_voxel',num2str(j),'.mat']);
    
    % remove excluded participants
    mvpa_data_seq = mvpa_data_seq(:,subsidx);
    mvpa_data_con = mvpa_data_con(:,subsidx);
    mvpa_data_congrap = mvpa_data_congrap(:,subsidx);
    
    % for each participant
    for i = 1:length(subs)
        % for each artefact location
        for a = 2:4    
            
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
            mvpa_collate_mean_seq(i,a-1) = mvpa_comparison_mean;

            % extract data for single-band control sequence
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
            mvpa_collate_mean_con(i,a-1) = mvpa_comparison_mean;

            % extract data for single-echo control sequence
            x = mvpa_data_congrap(a,i).data;
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
            mvpa_collate_mean_congrap(i,a-1) = mvpa_comparison_mean;
        end
    end
    
    % update violin plots
    
    % reformat results into a table for plotting
    tmp = zeros(17,6);
    tmp(:,logical([1 0 1 0 1 0])) = mvpa_collate_mean_seq;
    tmp(:,logical([0 0 0 1 0 0])) = mvpa_collate_mean_con(:,2);
    tmp(:,logical([0 1 0 0 0 1])) = mvpa_collate_mean_congrap(:,[1,3]);
    tmp = array2table(tmp,'VariableNames',roiname);
    % make a violin plot
    v = violinplot(tmp);
    % y-axis label
    ylabel('Difference in cosine distance');
    v(1,1).ViolinColor = colours(7,:);
    v(1,1).ViolinAlpha = 0.4;
    v(1,2).ViolinColor = colours(7,:);
    v(1,2).ViolinAlpha = 0.15;
    v(1,3).ViolinColor = colours(5,:);
    v(1,3).ViolinAlpha = 0.4;
    v(1,4).ViolinColor = colours(5,:);
    v(1,4).ViolinAlpha = 0.15;
    v(1,5).ViolinColor = colours(3,:);
    v(1,5).ViolinAlpha = 0.4;
    v(1,6).ViolinColor = colours(3,:);
    v(1,6).ViolinAlpha = 0.15;
    set(gca,'FontSize',16)
    ylim([-0.1,0.5])
    yticks(-0.1:0.1:0.5)
    xtickangle(gca,-45)
    set(gcf,'position',[10,10,600,600])
    saveas(gcf,['/group/mlr-lab/Saskia/7T/slice_leakage_figures/MEMB_peak',num2str(j),'_mvpa.png'])
    close(gcf)
    
    % update t-tests
    % A
    
    % A_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_congrap(:,1),mvpa_collate_mean_seq(:,1),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,1) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,1) = tmp4.tstat;
    
    % B
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_con(:,2),mvpa_collate_mean_seq(:,2),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,2) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,2) = tmp4.tstat;
    
    % B_g
    [tmp1,tmp2,tmp3,tmp4] = ttest(mvpa_collate_mean_congrap(:,3),mvpa_collate_mean_seq(:,3),'tail','left');
    % get p-values
    MEMB_mvpa_p_val(j,3) = tmp2;
    % also get t-values
    MEMB_mvpa_t_val(j,3) = tmp4.tstat;

end
