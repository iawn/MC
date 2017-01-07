function Cluster3DMC(directory,loadStart,numToLoad,rect,basename,parallelDepths)

% This code will motion correct and do pixel CI for all of the image stacks
% located in the specified directory
%clear all; close all; clc;
%file = '/global/scratch/mardinly/AM63_1a/20160303_trimOff/';
%loadStart=1;
%numToLoad=304;
%basename is a variable term if omitted will select all files in directory


%if dir is selected, make sure you provide a direcotry
if exist(directory)~=7;
    errordlg('Please provide directory')
end;

if nargin<5 || isempty(basename) || ~ischar(basename)
    basename ='';
    parallelDepths =0;
end

if nargin <6
    parallelDepths=0;
end


%% get all files in dir
fprintf('Identifying Files... \n'); 
k=dir(directory);k(1:2)=[];
i=1;
for n=1:(numel(k))
    fn=k(n).name;
    if regexp(fn,regexptranslate('wildcard',[basename '*.tif']))
        filenames{i}=fullfile(directory,  k(n).name);
        i=i+1;
    end;
end;
fprintf([ num2str(numel(filenames)) ' Files Detected... \n']);

%% Extract MetaData

s = imfinfo(filenames{1});
header = s(1).Software;
MD = parseSI5Header(header); 

%% set rectangle
try
    nDepths = MD.hFastZ.numFramesPerVolume;
catch
    nDepths = 1;
end


if ~exist('rect');
    rect=[1 1 511 511];
    rect=repmat(rect,[nDepths 1]);
    disp('no rectangle, using whole image')
elseif size(rect,2) ~=4
    rect=[1 1 511 511];
    rect=repmat(rect,[nDepths 1]);
    disp('rect did not have 4 elements, using whole image')
elseif size(rect,1) ~= nDepths
    rect=repmat(rect(1,:),[nDepths 1]);
    disp('rect only included for one depth, appying to all')
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




%% check and see if .align file exists already in dir
%  if not, run MC and possibly CI

nn = loadStart:(numToLoad+loadStart-1);
f=filenames(nn);

k=f{1}; k(length(k)-6:length(k))=[];
outputBase=[k num2str(loadStart) '-' num2str(numToLoad+loadStart-1)];
useTheseFiles=filenames(nn)';

% parallelDepths =0;

if parallelDepths
    
    
    nc = feature('numcores'); %number of cores
    availCores = nc -nDepths; %each depth takes a core to launch
    minCore = floor(availCores/nDepths);
    jobCores = ones(nDepths,1)*minCore;
    remCore = rem(availCores,nDepths);
    for i=1:remCore %distribute unused cores across jobs
        jobCores(i)=jobCores(i)+1;
    end
    
    
    for depth=1:nDepths;
        disp(['Launching Batch: ' num2str(depth)]);
        j(depth) = batch(['sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,' num2str(depth) ')'],...
            'AttachedFiles', 'sbxAlignOneDepth.m',...
            'Pool',jobCores(depth));
        
    end;
    
    %wait for jobs to end
    for i=1:numel(j)
        disp(['Waiting for job ' num2str(i) ' (others might still be running)']);
        wait(j(i));
        disp(['Diary for Depth ' num2str(i) ': ']);
        diary(j(i));
    end
    disp('jobs completed');
else
    for depth=1:nDepths;
        disp(['Starting Depth: ' num2str(depth)]);
        sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,depth)
    end
    disp('job completed');
end


end


