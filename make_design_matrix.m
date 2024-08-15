names = cell(1,2);
onsets = cell(1,2);
durations = cell(1,2);

names{1} = 'Semantic';
onsets{1} = [16    48    96   128   176   208   256   288   336   368   416   448];
durations{1} = 16;

names{2} = 'Control';
onsets{2} = [32    64   112   144   192   224   272   304   352   384   432   464];
durations{2} = 16;

save('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/design_matrix.mat','names','onsets','durations');