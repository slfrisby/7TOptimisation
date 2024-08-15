% summarises accuracy and reaction time and enables identification of
% participants who do not perform as expected.

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

% ACCURACY

% calculate percentage accuracy for each run for each participant
for i=1:5
    percentage_accuracy(i,:)=((sum(behav_reordered{i}.Acc,1)/96)*100);
end

% (note that accuracy is not normally distributed)
% figure
% histogram(percentage_accuracy)
% therefore, use median rather than mean.

% this is a repeated-measures design, so the strategy is to calculate a
% summary statistic PER PARTICIPANT and inspect the distribution of that
% summary statistic across participants.
median_accuracy = nanmedian(percentage_accuracy);

% plot
figure
boxplot(median_accuracy)
title('Accuracy')
% there are no visible outliers.

% REACTION TIME

% we only want to calculate reaction time for the trials on which the
% participant answered correctly.
for i = 1:5
    tmp = behav_reordered{i}.RT;
    mask = behav_reordered{i}.Acc;
    tmp(mask==0)=NaN;
    reaction_time(i,:) = nanmean(tmp);
end

% reaction time is normally distributed, but use median for consistency
% figure
% histogram(reaction_time)

median_rt = nanmedian(reaction_time);

% plot
figure
boxplot(median_rt)
title('Reaction Time')
% there are no visible outliers.





