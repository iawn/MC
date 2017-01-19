function sbxalignmastermulti_3D(ImageFile,computeci,outputname,hi,wi,MD,depth,useRed)

szz(1)=numel(hi);
szz(2)=numel(wi);

nfil=length(ImageFile);
nDepth = MD.hFastZ.numFramesPerVolume;
numVolumes =MD.hFastZ.numVolumes;
nFramesTot = nDepth * numVolumes * 2; %total number of frames in an acq including red

%Computing first order stats
lTime = tic;
fprintf('Getting first-order stats\n');
IF=1:nfil; %image index
if useRed
    IFr=2:2:2*MD.hFastZ.numVolumes * MD.hFastZ.numFramesPerVolume; % get red frame indx
else
    IFr=1:2:2*MD.hFastZ.numVolumes * MD.hFastZ.numFramesPerVolume; % get green frame indx
end
IFr=IFr(depth:nDepth:end); %get just this depth

theMatrix=combvec(IFr,IF);
theMatrix=theMatrix';

ms = 0;
vs = 0;
%s=size(theStack);
X = [ones(length(theMatrix),1),linspace(-1,1,length(theMatrix))'];
X = bsxfun(@times,X,1./sqrt(sum(X.^2)));

%Find appropriate chunking volume to process some frames at a time;
maxChunk = 500;
minChunk = 50; %will overwrite if can't find a better min

if numVolumes <=maxChunk;
    Chunk = numVolumes;
else

M=2;
p = perms(factor(numVolumes)) ; 
y = unique(sort([p(:,1:M-1) prod(p(:,M:end),2)],2),'rows') ; %list of viable chunks ie 2 chunks of 300 frames

Chunk = y(find(y(:,2)<=maxChunk & y(:,2)>=minChunk,1),2); %to load 50 to 500 frames at a time;
while isempty(Chunk) && M<5
    M=M+1; %if it is a large number of frames in each acq this will help
    y = unique(sort([p(:,1:M-1) prod(p(:,M:end),2)],2),'rows') ; %list of viable chunks ie 2 chunks of 300 frames
    Chunk = y(find(y(:,M)<=maxChunk & y(:,M)>=minChunk,1),M); %to load 50 to 500 frames at a time;
end

if isempty(Chunk) %if its still empty
    fact = factor(numVolumes);
    Chunk = fact(end);
    fprintf(['could not find a good frame chunking number going at ' num2str(Chunk) '\n']);
end
end

numChunks = MD.hFastZ.numVolumes / Chunk;


%ImageFile{theMatrix(jj,2)}
effChunks = Chunk * 2 * nDepth; %effective chunks =frames *2 colors * n depths
rc1 =0;
for i=1:numChunks*nfil

    [~,~,z] = bigread3(ImageFile{theMatrix((i-1)*Chunk+1,2)},...
        rem(1+(i-1)*effChunks,nFramesTot),...
        effChunks);%load2P(ImageFile,'Frames',jj,'Double');
   
        z = z(:,:,1:2:end);%select color (should be ever other from the start color regardless if thats green or red 
    
z = z(:,:,depth:nDepth:end); %just this depth

parfor jj = 1:Chunk

    z2 = double(z(:,:,jj));

    z2 = z2(hi,wi);

    %z=theStack(:,:,jj);

    ms = ms + z2(:)*X(jj+(i-1)*Chunk,:);

    vs = vs + z2(:).^2;
    rc1=rc1+1;
end

end

% %%%%temp debug stuff
% ms2 = 0;
% vs2 = 0;
% rc2 =0;
% parfor ii = 1:length(theMatrix)
%     zx = bigread3(ImageFile{theMatrix(ii,2)},...
%         theMatrix(ii,1),...
%         1);
%      z2 = double(zx);
% 
%     z2 = z2(hi,wi);
% 
%     %z=theStack(:,:,jj);
% 
%     ms2 = ms2 + z2(:)*X(ii,:);
% 
%     vs2 = vs2 + z2(:).^2;
%     rc2=rc2+1;
% end


fprintf(['First Order Stats done. Time ' num2str(toc(lTime)) 's \n'])
k=sqrt(1/MD.hFastZ.numVolumes *(vs - sum(ms.^2,2)));
s = reshape(k,[szz]);


try
    thestd = medfilt2(s,[31,31],'symmetric');
catch
    thestd = medfilt2(real(s),[31,31],'symmetric');
end



gl = X(:,2);

l  = reshape(ms(:,2),szz);



%%

%for imgFileN=1:nImgFile;
%fprintf(['Working on Stack #' num2str(imgFileN) '\n']');
fprintf('Alignment first pass\n');

%ms=msi;
%vs=vsi;
%gl=gli;
%l=li;
%s=si;
%thestd=thestdi;

% %%%debug stuff
% l2  = reshape(ms2(:,2),szz);
% k2=sqrt(1/MD.hFastZ.numVolumes *(vs2 - sum(ms2.^2,2)));
% s2 = reshape(k2,[szz]);
% thestd2 = medfilt2(s2,[31,31],'symmetric');

sapmTime  = tic;
[m,~,T] = sbxalignparmulti6(ImageFile,thestd,gl,l,theMatrix,hi,wi); %Takes about 2.5 minutes
% %debug stuff
% [m2,~,T2] = sbxalignparmulti6(ImageFile,thestd2,gl,l2,theMatrix,hi,wi); %Takes about 2.5 minutes
fprintf(['First alignment ' num2str(toc(sapmTime)) 's\n']);

rgx = (1:size(m,2))+45;

rgy = 32 + (1:size(m,1));

T0 = T;



for nn = 1:10
    passTime = tic;
    fprintf('Refining alignment... pass %d\n',nn);
    
    [m,~,T] = sbxaligniterativemulti6(ImageFile,m,rgy,rgx,thestd(rgy,rgx),gl,l,theMatrix,hi,wi);
    
    dT = sqrt(mean(sum((T0-T).^2,2)));
    
    T0 = T;
    
    fprintf('delta: %.3f\n',dT);
    fprintf(['This pass ' num2str(toc(passTime)) 's\n']); 
    
    if dT < .25
        
        break;
        
    end
    
    
    
end



fprintf('Getting aligned first-order stats\n');
alignFOStime = tic;

ms = 0;

vs = 0;



m2 = 0;

v2 = 0;



X = [ones(length(theMatrix),1),linspace(-1,1,length(theMatrix))'];

X = bsxfun(@times,X,1./sqrt(sum(X.^2)));



g = exp(-(-5:5).^2/2/1.6^2);




for i=1:numChunks*nfil;

     z1 = bigread3(ImageFile{theMatrix((i-1)*Chunk+1,2)},...
        rem(1+(i-1)*effChunks,nFramesTot) + theMatrix(1,1)-1,... %include offset for red channel or skipped early frames **IMPORTANT** assumes every file has the same offset
        effChunks);%load2P(ImageFile,'Frames',jj,'Double');

%z1 = z1(:,:,1:2:end);%just this color (doesn't matter green or red)
z1 = z1(:,:,depth:nDepth:end); %just this depth

parfor jj = 1:Chunk;

        z = double(z1(:,:,jj));
         z = z(hi,wi);
    
   %z = single(load2P(ImageFile,'Frames',jj));
   %z=z(:,:,1,1);
   %z=theStack(:,:,jj);
   z = z./thestd;
    
    z = circshift(z,T(jj+(i-1)*Chunk,:));
    
    z = double(z);
    
    
    
    ms = ms + z(:)*X(jj+(i-1)*Chunk,:);
    
    vs = vs + z(:).^2;
    
    
    
    z = conv2(g,g,z,'same');
    
    m2 = m2 + z(:)*X(jj+(i-1)*Chunk,:);
    
    v2 = v2 + z(:).^2;
    
end
end

fprintf(['Aligned First Order Stat Time: ' num2str(toc(alignFOStime)) 's\n']);


ss = sqrt(1/length(theMatrix)*(vs - sum(ms.^2,2)));

m = reshape(ms(:,1)*X(1,1),size(l));

v = reshape(ss.^2,size(l));





[Q,~] = qr(ms,0);

[Q2,~] = qr(m2,0);

s2 = reshape(sqrt(1/length(theMatrix)*(v2 - sum(m2.^2,2))),size(l));

m2 = reshape(m2(:,1)*X(1,1),size(l));






k = 0;

fprintf('Computing simple stats... pass %d\n',2);
simpleStatsTime = tic;


for i=1:numChunks*nfil;

    z1 = bigread3(ImageFile{theMatrix((i-1)*Chunk+1,2)},...
        rem(1+(i-1)*effChunks,nFramesTot)+ theMatrix(1,1)-1,... %include offset for red channel or skipped early frames **IMPORTANT** assumes every file has the same offset
        effChunks);%load2P(ImageFile,'Frames',jj,'Double');

%z1 = z1(:,:,1:2:end);%just this color (doesn't matter green or red)
z1 = z1(:,:,depth:nDepth:end); %just this depth

parfor jj = 1:Chunk;

        z = double(z1(:,:,jj));
         z = z(hi,wi);
   % z = load2P(ImageFile,'Frames',jj,'Double');
    %z = z(:,:,1,1);
  %  z=theStack(:,:,jj);

    z = z - gl(jj)*l;
    
    
    z = circshift(z,T(jj+(i-1)*Chunk,:));
    
    z = z./thestd;
    
    
    
    z = conv2(g,g,z,'same');
    
    z = reshape(z(:) - (Q2*(Q2'*z(:))),size(z));
    
    k  =  k + (z./s2).^4;
    
end
end

fprintf(['Simple Stats Time: ' num2str(toc(simpleStatsTime)) 's\n']);


sm = s2./m2;

k = k/length(theMatrix) - 3;


%[p,f,~] = fileparts(ImageFile{1});
%fname = fullfile(p,f);

try
%    save([fname '_align.mat'],'m','v','thestd','sm','k','T','Q','-mat');

    save(outputname,'m','v','thestd','sm','k','T','Q','-mat','-append');
    %save([fname '.alignDebug'],'msi','vsi','thestdi','gli','li','si','-mat','-append');

catch
    
    save(outputname,'m','v','thestd','sm','k','T','Q','-mat','-v7.3');
   % save([fname '.alignDebug'],'msi','vsi','thestdi','gli','li','si','-mat','-v7.3');

end


 
if computeci
fprintf(['Starting CI\n']');
ciTime = tic;
sbxcomputeci3D(ImageFile,outputname,hi,wi,MD); %Takes about 10 minutes, eats up a ton of RAM
fprintf(['CI time: ' num2str(toc(ciTime)) 's\n']);

end 

%end;
