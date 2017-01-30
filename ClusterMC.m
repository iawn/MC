function ClusterMC(file,basename,loadStart,numToLoad,rect)

% This code will motion correct and do pixel CI for all of the image stacks
% located in the specified directory

%% Arg Checks


%if dir is selected, make sure you provide a direcotry
if exist(file)~=7;
    errordlg('Please provide directory')
end;


%% get all files in dir
fprintf('Identifying Files... \n'); 
k=dir(file);k(1:2)=[];
i=1;
for n=1:(numel(k))
    fn=k(n).name;
    if regexp(fn,regexptranslate('wildcard',[basename '*.tif']))
        filenames{i}=fullfile(file,  k(n).name);
        i=i+1;
    end;
end;
fprintf([ num2str(numel(filenames)) ' Files Detected... \n']);

%% Determine Rectangle

if ~exist('rect');
    rect(1)=1;
    rect(2)=1;
    rect(3)=511;
    rect(4)=511;
    fprintf('no rectangle, using whole image\n')
elseif numel(rect)~=4;
    rect(1)=1;
    rect(2)=1;
    rect(3)=511;
    rect(4)=511;
    fprintf('rect did not have 4 elements, using whole image\n')
else
    fprintf('Rect Found...\n')
end;

wi=rect(1):rect(1)+rect(3);
hi=rect(2):rect(2)+rect(4);

rectMC = rect;
%% check and see if .align file exists already in dir
%  if not, run MC and possibly CI

nn = loadStart:(numToLoad+loadStart-1);
f=filenames(nn);

k=f{1}; k(length(k)-6:length(k))=[];
outputBase=[k num2str(loadStart) '-' num2str(numToLoad+loadStart-1)];

if ~exist([outputBase '.align']);
    useTheseFiles=filenames(nn)';
    fprintf('Beginning MC\n')
    sbxalignmastermulti6(useTheseFiles,1,[outputBase '.align'],hi,wi);
else
    fprintf('Align File Already Exists! \n');
end;

fprintf('Saving... \n')
save([outputBase '.align'],'outputBase','-append');
save([outputBase '.align'],'file','-append');
save([outputBase '.align'],'loadStart','-append');
save([outputBase '.align'],'numToLoad','-append');
save([outputBase '.align'],'filenames','-append');
save([outputBase '.align'],'rectMC','-append');



