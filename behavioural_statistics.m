% test whether there are any differences in accuracy or reaction time
% between sequences or between semantic and control tasks.

cd('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/supplementary');

% loop over participants
for s=1:20

% data is copied from EPrime output into an Excel spreadsheet, one sheet
% per participant. Each sheet is labelled with the participant number:
subcode{s}=num2str(s,'%03.f');
% read the correct sheet
x=xlsread('accuracy_and_rt.xlsx',subcode{s});

% read accuray and reaction time from that sheet and store them in a cell
% array. The cell array contains 1 cell array for each participant, which
% contains 1 struct per run, which contains 2 fields - reaction time and
% accuracy.
%run1
behav{s}{1}.RT=x(1:96,1);
behav{s}{1}.Acc=x(1:96,4);
%run2
behav{s}{2}.RT=x(1:96,8);
behav{s}{2}.Acc=x(1:96,11);
%run3
behav{s}{3}.RT=x(1:96,15);
behav{s}{3}.Acc=x(1:96,18);
%run4
behav{s}{4}.RT=x(1:96,22);
behav{s}{4}.Acc=x(1:96,25);
%run5
behav{s}{5}.RT=x(1:96,29);
behav{s}{5}.Acc=x(1:96,32);

% Reorder the data into another cell array containing 1 struct per run,
% containing 2 fields (reaction time and accuracy), containing accuracy and
% reaction time collated into one matrix.
behav_reordered{1}.RT(:,s)=x(1:96,1);
behav_reordered{2}.RT(:,s)=x(1:96,8);
behav_reordered{3}.RT(:,s)=x(1:96,15);
behav_reordered{4}.RT(:,s)=x(1:96,22);
behav_reordered{5}.RT(:,s)=x(1:96,29);

behav_reordered{1}.Acc(:,s)=x(1:96,4);
behav_reordered{2}.Acc(:,s)=x(1:96,11);
behav_reordered{3}.Acc(:,s)=x(1:96,18);
behav_reordered{4}.Acc(:,s)=x(1:96,25);
behav_reordered{5}.Acc(:,s)=x(1:96,32);

end

% for every run, the pattern of stimuli was 4 semantic, 4 control, 4
% semantic, 4 control etc. Make a filter to filter between tasks
taskmask = [0;0;0;0;1;1;1;1];
taskmask = repmat(taskmask,12,20);

% each participant underwent each imaging sequence in a counterbalanced
% order. 
counterbalancing = readcell('counterbalancing.xlsx');
counterbalancing = counterbalancing(2:21,7:11);

% unscramble behav_reordered so that all SESB values are in column 1, all SEMB in column 2,
% all MESB in column 3, all MEMB in column 4, all pTx in column 5
for s=1:20
    for i=1:5
        switch counterbalancing{s,i}
            case 'standard'
                behav_sequence{1}.RT(:,s) = behav_reordered{i}.RT(:,s);
                behav_sequence{1}.Acc(:,s) = behav_reordered{i}.Acc(:,s);
            case 'multiband'
                behav_sequence{2}.RT(:,s) = behav_reordered{i}.RT(:,s);
                behav_sequence{2}.Acc(:,s) = behav_reordered{i}.Acc(:,s);
            case 'multiecho'
                behav_sequence{3}.RT(:,s) = behav_reordered{i}.RT(:,s);
                behav_sequence{3}.Acc(:,s) = behav_reordered{i}.Acc(:,s);
            case 'MBME'
                behav_sequence{4}.RT(:,s) = behav_reordered{i}.RT(:,s);
                behav_sequence{4}.Acc(:,s) = behav_reordered{i}.Acc(:,s);
            case 'pTx'
                behav_sequence{5}.RT(:,s) = behav_reordered{i}.RT(:,s);
                behav_sequence{5}.Acc(:,s) = behav_reordered{i}.Acc(:,s);
            otherwise % to catch typos!
                disp(['Suspected typo! s = ',num2str(s), ', i = ',num2str(i)]);
        end
    end
end

% ACCURACY 

% loop over runs
for i=1:5
    % get accuracy for that run for all participants
    tmp = behav_sequence{i}.Acc;
    % semantic task first - set all scores from the control task to NaN
    tmp(taskmask==1) = NaN;
    accuracy{1}{i} = ((nansum(tmp,1)/48)*100);
    % then control task
    tmp = behav_sequence{i}.Acc;
    tmp(taskmask==0) = NaN;
    accuracy{2}{i} = ((nansum(tmp,1)/48)*100);
end

% construct summary table of median and IQR
for i = 1:5
    for task = 1:2
        median_accuracy(task,i) = nanmedian(accuracy{task}{i});
        iqr_accuracy(task,i) = iqr(accuracy{task}{i});
    end
end

% Friedman's test for differences between sequences (columns)
friedman(median_accuracy)
% t-test for differences between tasks (rows)
[p,h,stats] = signrank(median_accuracy(1,:),median_accuracy(2,:))

% REACTION TIME

for i=1:5
    % get reaction time for that run for all participants
    tmp = behav_sequence{i}.RT;
    % mask for correct trials only
    tmp(behav_sequence{i}.Acc==0) = NaN;
    % semantic task first - set all scores from the control task to NaN
    tmp(taskmask==1) = NaN;
    reaction_time{1}{i} = nanmean(tmp);
    % then control task
    tmp = behav_sequence{i}.RT;
    tmp(behav_sequence{i}.Acc==0) = NaN;
    tmp(taskmask==0) = NaN;
    reaction_time{2}{i} = nanmean(tmp);  
end

% construct summary table of median and IQR
for i = 1:5
    for task = 1:2
        median_reaction_time(task,i) = nanmedian(reaction_time{task}{i});
        iqr_reaction_time(task,i) = iqr(reaction_time{task}{i});
    end
end

% Friedman's test for differences between sequences (columns)
friedman(median_reaction_time)
% t-test for differences between tasks (rows)
[p,h,stats] = signrank(median_reaction_time(1,:),median_reaction_time(2,:))
