addpath('/group/mlr-lab/AH/Projects/spm12/');
addpath('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts')

%load atlas .nii
x=spm_vol('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/Comb_atlas.nii');
x_data=spm_read_vols(x);
% load atlas labels
atlas_label = readtable('/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/scripts/Comb_atlas.xlsx','ReadVariableNames',false);

% initialise output - mni coordinates as matrix x,y,z columns
spm_table=[];
% initialise output - labels
spm_label={};

% fill in table of mni coordinates
for i=1:size(TabDat.dat,1)
    tmp=TabDat.dat{i,12}';
    spm_table(i,:)=tmp;
end

% convert mni coordinates into matrix space. Find the location of this
% point within the atlas. The numeric value of the atlas at that point
% gives the anatomical label. 
for i=1:size(spm_table,1)
    %convert to matrix space
    t=mni2cor(spm_table(i,:),x.mat);
    %extract coordinate value in matrix space
    if isempty(find([t(1),t(2),t(3)]<1)) % unless the coordinates lie outside the atlas (in which case you will have to open the atlas in mricron, find the coordinate, and find the nearest region)
        out(i,:)=x_data(t(1),t(2),t(3));
        % get label, including hemisphere label
        if out(i,:)==0 % if the coordinate is outside any regions in the atlas
            spm_label(i,:)={''}; % leave the row empty - you will then have to open the atlas in mricron, find the coordinate, and find the nearest region
        else
            spm_label(i,:)={[char(atlas_label{out(i,:),3}),' ',lower(char(atlas_label{out(i,:),2}))]};
        end
    end
end


