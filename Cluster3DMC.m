function Cluster3DMC(directory,loadStart,numToLoad,rect,basename,useRed,useDepths,parallelDepths)
% This code will motion correct and do pixel CI for all of the image stacks
% located in the specified directory
%
%   Directory       - is the path to all of the files to be motion corrected
%   loadStart       - Which file to begin the loading on (set to 1 for all
%                     files)
%   numToLoad       - Number of files to Motion Correct
%   rect            - (optional) rectangle to exclude edges from motion
%                     correction. defaults to [1 1 511 511]
%   basename        - (optional) basename that all Motion Corrected Files must
%                     match
%   useRed          - (optional) align using red channel.
%   useDepths       - (optional) array of depths to be processed this run.
%                     defaults to all.
%   parallelDepths  - (optional) different form of parallelization that
%                     seperates each depth. In practice doesn't save much 
%                     time as the load steps are the slow ones.



%if dir is selected, make sure you provide a direcotry
if exist(directory)~=7
    errordlg('Please provide directory')
end;

if nargin<4
    rect=[1 1 511 511];
    basename ='';
    parallelDepths =0;
    useRed =0;
    useDepths=0;
elseif nargin<5 
    basename ='';
    parallelDepths =0;
    useRed =0;
    useDepths=0;
elseif nargin<6
    useRed =0;
    parallelDepths=0;
    useDepths=0;
elseif nargin <7
    useDepths=0;
    parallelDepths=0;
elseif nargin <8
    parallelDepths=0;
end

if isempty(basename) || ~ischar(basename)
  basename='';
end

verboseFlag = 0;

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
    %if useDepths~=0;
    %   nDepths=numel(useDepths);
    %end
catch
    nDepths = 1;
end

if useDepths == 0
    useDepths = 1:nDepths;
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
    
    jobCores=max(jobCores,1);
    
    disp([num2str(availCores) ' cores available. divided as: ' num2str(jobCores')]);
    
    for depth=useDepths
        disp(['Launching Batch: ' num2str(depth)]);
        j(depth) = batch(['sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,' num2str(depth) ',useRed)'],...
            'AttachedFiles', 'sbxAlignOneDepth.m',...
            'Pool',jobCores(depth));
    end
    
    
    
    %wait for jobs to end
    for i=1:numel(j)
        jTime = tic;
        disp(['Waiting for Job ' num2str(i) ' (others might still be running)']);
        wait(j(i));
        disp(['Depth ' num2str(i) ' took ' num2str(toc(jTime)) 's additional']);
        if verboseFlag
            disp(['Diary for Depth ' num2str(i) ': ']);
            diary(j(i));
        end

    end
    delete(j);
    clear j;

    for depth=useDepths
        OPB=[outputBase '_depth_' num2str(i)];
        fprintf('Saving Files\n');
        
        save([OPB '.align'],'outputBase','-append');
        save([OPB '.align'],'useTheseFiles','-append');
        save([OPB '.align'],'loadStart','-append');
        save([OPB '.align'],'numToLoad','-append');
        save([OPB '.align'],'hi','wi','-append');
        try
            save([OPB '.align'],'rect','-append');
        end
        
    end;
    
    disp('jobs completed');
else
    poolTime = tic;
    nc = feature('numcores'); %number of cores
    disp([num2str(nc) ' cores available.']);
    try
        parpool(nc);
    catch
        fprintf('Could not launch the desired number of cores switching to default\n')
        parpool
    end
    fprintf(['Launching parpool took ' num2str(toc(poolTime)) 's\n']);
    
    for depth=useDepths
        %thisDepth=useDepths(depth);
        dTime=tic;
        disp(['Starting Depth: ' num2str(depth)]);
        sbxAlignOneDepth(outputBase,useTheseFiles,hi,wi,MD,depth,useRed)
        disp(['Depth ' num2str(depth) ' took ' num2str(toc(dTime)) 's']);
        
        OPB=[outputBase '_depth_' num2str(depth)];
        fprintf('Saving Files\n');
        save([OPB '.align'],'outputBase','-append');
        save([OPB '.align'],'useTheseFiles','-append');
        save([OPB '.align'],'loadStart','-append');
        save([OPB '.align'],'numToLoad','-append');
        save([OPB '.align'],'hi','wi','-append');
        try
            save([OPB '.align'],'rect','-append');
        end
    end
    
    disp('job completed');
end


end


