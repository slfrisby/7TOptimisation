% calculates, for each participant and each run, the percentage of volumes
% in that run for which motion was greater than 2mm translation and/or 1
% degree rotation.

% set root folder
root='/imaging/projects/cbu/wbic-p00567-7Tmultiecho/main/work';
cd(root)
% get all folders starting with "sub-0"
folders=dir('sub-*');

%run labels
cond=[{'MESB'},{'MEMB'},{'SESB'},{'SEMB'},{'ptx8ms'}];

%initialise output
outlier=zeros(length(folders),length(cond));

%loop subjects
for s=1:length(folders)
    %loop runs
    for i=1:length(cond)
        % no pTx file for participant 4, so skip over this one
        if (s==4) && (i==5)
            continue
        end
        %cd to folder
        cd([root,'/',folders(s).name,'/',cond{i}]);
        %check if motion file exists
        if exist('ratmpdata.1D')
            %load file if does
            x=load('ratmpdata.1D');
            
            %binarise any absolute rotation values greater than 1 degree
            rot=double((abs(x(:,2:4)))>1);
            %binarise any absolute translation values greater than 2 mm
            tran=double((abs(x(:,5:7)))>2);
            %combine both
            motion=[rot,tran];
            %get sum for any volume that is an outlier
            motion=sum(motion,2);
            %calculate % of outliers based on total run length
            %save in outlier subject x run matrix
            outlier(s,i)=(sum(motion>0.5)/size(x,1))*100;
        end

    end

end

%print variable to screen to see outliers and plot figure to visualise
outlier
figure;
subplot(2,1,1);
imagesc(outlier',[0 100]);
subplot(2,1,2);
boxplot(outlier');
