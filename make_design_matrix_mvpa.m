names = cell(1,24);
onsets = cell(1,24);
durations = cell(1,24);

names{1} = 'Semantic1';
onsets{1} = [16];
durations{1} = [16];

names{2} = 'Semantic2';
onsets{2} = [48];
durations{2} = [16];

names{3} = 'Semantic3';
onsets{3} = [96];
durations{3} = [16];

names{4} = 'Semantic4';
onsets{4} = [128];
durations{4} = [16];

names{5} = 'Semantic5';
onsets{5} = [176];
durations{5} = [16];

names{6} = 'Semantic6';
onsets{6} = [208];
durations{6} = [16];

names{7} = 'Semantic7';
onsets{7} = [256];
durations{7} = [16];

names{8} = 'Semantic8';
onsets{8} = [288];
durations{8} = [16];

names{9} = 'Semantic9';
onsets{9} = [336];
durations{9} = [16];

names{10} = 'Semantic10';
onsets{10} = [368];
durations{10} = [16];

names{11} = 'Semantic11';
onsets{11} = [416];
durations{11} = [16];

names{12} = 'Semantic12';
onsets{12} = [448];
durations{12} = [16];

names{13} = 'Control1';
onsets{13} = [32];
durations{13} = [16];

names{14} = 'Control2';
onsets{14} = [64];
durations{14} = [16];

names{15} = 'Control3';
onsets{15} = [112];
durations{15} = [16];

names{16} = 'Control4';
onsets{16} = [144];
durations{16} = [16];

names{17} = 'Control5';
onsets{17} = [192];
durations{17} = [16];

names{18} = 'Control6';
onsets{18} = [224];
durations{18} = [16];

names{19} = 'Control7';
onsets{19} = [272];
durations{19} = [16];

names{20} = 'Control8';
onsets{20} = [304];
durations{20} = [16];

names{21} = 'Control9';
onsets{21} = [352];
durations{21} = [16];

names{22} = 'Control10';
onsets{22} = [384];
durations{22} = [16];

names{23} = 'Control11';
onsets{23} = [432];
durations{23} = [16];

names{24} = 'Control12';
onsets{24} = [464];
durations{24} = [16];

save('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/design_matrix_mvpa.mat','names','onsets','durations');