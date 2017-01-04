function Cluster3DMC(file,loadStart,numToLoad,rect,nDepths);

% This code will motion correct and do pixel CI for all of the image stacks
% located in the specified directory
%clear all; close all; clc;
%file = '/global/scratch/mardinly/AM63_1a/20160303_trimOff/';
%loadStart=1;
%numToLoad=304;

if ~exist('rect');
rect=[1 1 511 511];
rect=repmat(rect,[nDepths 1]);
disp('no rectangle, using whole image')
elseif numel(rect(:,1))~=nDepths;
rect=[1 1 511 511];
rect=repmat(rect,[nDepths 1]);
disp('rect did not have 4 elements, using whole image')
end;

for n=1:nDepths
wiTemp=rect(n,1):rect(n,1)+rect(n,3);
hiTemp=rect(n,2):rect(n,2)+rect(n,4);

%check Range Errors
wiTemp(wiTemp>=511)=[];
wiTemp(wiTemp<=1)=[];
hiTemp(hiTemp>=511)=[];
hiTemp(hiTemp<=1)=[];

wi{n}=wiTemp;
hi{n}=hiTemp;

end



%if dir is selected, make sure you provide a direcotry
if exist(file)~=7;
    errordlg('Please provide directory')
end;




%% get all files in dir
k=dir(file);k(1:2)=[];
i=1;
for n=1:(numel(k))
    fn=k(n).name;
    if strcmp(fn(length(fn)-2:length(fn)),'tif')
        filenames{i}=[file  k(n).name];
        i=i+1;
    end;
end;

%% Extract MetaData

s = imfinfo(filenames{1});
header = s(1).Software;
MD = parseSI5Header(header); 


%% check and see if .align file exists already in dir
%  if not, run MC and possibly CI

nn = loadStart:(numToLoad+loadStart-1);
f=filenames(nn);

k=f{1}; k(length(k)-6:length(k))=[];
outputBase=[k num2str(loadStart) '-' num2str(numToLoad+loadStart-1)];
useTheseFiles=filenames(nn)';

for depth=1:nDepths;

OPB=[outputBase '_depth_' num2str(depth)];
disp('Beginning MC')
sbxalignmastermulti_3D(useTheseFiles,1,[OPB '.align'],hi{depth},wi{depth},MD,depth);

save([OPB '.align'],'outputBase','-append');
save([OPB '.align'],'file','-append');
save([OPB '.align'],'loadStart','-append');
save([OPB '.align'],'numToLoad','-append');
save([OPB '.align'],'hi','wi','-append');
try
save([OPB '.align'],'rect','-append');
end

%sbxcomputeci(useTheseFiles,[outputBase '.align'],hi{depth},wi{depth}); %Takes about 10 minutes, eats up a ton of RAM

end;
