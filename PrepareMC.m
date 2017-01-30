function PrepareMC(f,fn,saveloc)

%this function loads a single volume acq file from SI5 3D and returns a
%rectangle and metadata

if nargin<2;
    disp('specify file directory and savename')
    return;
end

if ~ischar(fn);
    disp('please specify file name as string')
    return;
end


Ims = ScanImageTiffReader(f).data();
meta = ScanImageTiffReader(f).metadata();
SI = parseSI5Header(meta);

nDepths=numel(SI.hStackManager.zs);


imG=Ims(:,:,1:2:end);
imR=Ims(:,:,2:2:end);

for j = 1:nDepths;
    close all    
    disp(['Draw Box around non-artifact area for depth ' num2str(j)])
    imagesc(mean(imR(:,:,j:nDepths:end),3));
    r = getrect;
    rect(j,:)=round(r);
end

if exist('saveloc')
save([saveloc fn],'SI','rect','nDepths');
else
save([fn],'SI','rect','nDepths');
end


