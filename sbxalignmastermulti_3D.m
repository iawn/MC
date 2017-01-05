function sbxalignmastermulti_3D(ImageFile,computeci,outputname,hi,wi,MD,depth)

szz(1)=numel(hi);
szz(2)=numel(wi);

nfil=length(ImageFile);
nDepth = MD.hFastZ.numFramesPerVolume;
numVolumes =MD.hFastZ.numVolumes;
nFramesTot = nDepth * numVolumes * 2; %total number of frames in an acq including red

%Computing first order stats

fprintf('Getting first-order stats\n');
IF=1:nfil; %image index
IFr=1:2:2*MD.hFastZ.numVolumes * MD.hFastZ.numFramesPerVolume; % get green frame indx
IFr=IFr(depth:nDepth:end); %get just this depth

theMatrix=combvec(IFr,IF);
theMatrix=theMatrix';

ms = 0;
vs = 0;
%s=size(theStack);
X = [ones(length(theMatrix),1),linspace(-1,1,length(theMatrix))'];
X = bsxfun(@times,X,1./sqrt(sum(X.^2)));

%find appropriate chunking volume to process some frames at a time;


maxChunk = 300;
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
    disp(['could not find a good frame chunking number going at ' num2str(Chunk)]);
end
end

numChunks = MD.hFastZ.numVolumes / Chunk;


%ImageFile{theMatrix(jj,2)}
effChunks = Chunk * 2 * nDepth; %effective chunks =frames *2 colors * n depths

for i=1:numChunks*nfil;

    z = bigread3(ImageFile{theMatrix((i-1)*Chunk+1,2)},...
        rem(1+(i-1)*effChunks,nFramesTot),...
        effChunks);%load2P(ImageFile,'Frames',jj,'Double');

%z = z(:,:,1:2:end);%green
z = z(:,:,depth:nDepth:end); %just this depth

for jj = 1:Chunk;

    z2 = double(z(:,:,jj));

    z2 = z2(hi,wi);

    %z=theStack(:,:,jj);

    ms = ms + z2(:)*X(jj,:);

    vs = vs + z2(:).^2;

end

end

disp('loaded')
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



tic
[m,~,T] = sbxalignparmulti6(ImageFile,thestd,gl,l,theMatrix,hi,wi); %Takes about 2.5 minutes
toc


rgx = (1:size(m,2))+45;

rgy = 32 + (1:size(m,1));

T0 = T;


tic
for nn = 1:10
    
    fprintf('Refining alignment... pass %d\n',nn);
    
    [m,~,T] = sbxaligniterativemulti6(ImageFile,m,rgy,rgx,thestd(rgy,rgx),gl,l,theMatrix,hi,wi);
    
    dT = sqrt(mean(sum((T0-T).^2,2)));
    
    T0 = T;
    
    if dT < .25
        
        break;
        
    end
    
    fprintf('delta: %.3f\n',dT);
    
end
toc


fprintf('Getting aligned first-order stats\n');



ms = 0;

vs = 0;



m2 = 0;

v2 = 0;



X = [ones(length(theMatrix),1),linspace(-1,1,length(theMatrix))'];

X = bsxfun(@times,X,1./sqrt(sum(X.^2)));



g = exp(-(-5:5).^2/2/1.6^2);

tic;

parfor jj = 1:length(theMatrix)
    
    
    
         z = bigread3(ImageFile{theMatrix(jj,2)},theMatrix(jj,1),1);%load2P(ImageFile,'Frames',jj,'Double');
         z = double(z);
         z = z(hi,wi);
    
   %z = single(load2P(ImageFile,'Frames',jj));
   %z=z(:,:,1,1);
   %z=theStack(:,:,jj);
   z = z./thestd;
    
    z = circshift(z,T(jj,:));
    
    z = double(z);
    
    
    
    ms = ms + z(:)*X(jj,:);
    
    vs = vs + z(:).^2;
    
    
    
    z = conv2(g,g,z,'same');
    
    m2 = m2 + z(:)*X(jj,:);
    
    v2 = v2 + z(:).^2;
    
end

toc;



ss = sqrt(1/length(theMatrix)*(vs - sum(ms.^2,2)));

m = reshape(ms(:,1)*X(1,1),size(l));

v = reshape(ss.^2,size(l));





[Q,~] = qr(ms,0);

[Q2,~] = qr(m2,0);

s2 = reshape(sqrt(1/length(theMatrix)*(v2 - sum(m2.^2,2))),size(l));

m2 = reshape(m2(:,1)*X(1,1),size(l));






k = 0;

fprintf('Computing simple stats... pass %d\n',2);

parfor jj = 1:length(theMatrix)
    
    
         z = bigread3(ImageFile{theMatrix(jj,2)},theMatrix(jj,1),1);%load2P(ImageFile,'Frames',jj,'Double');
         z = double(z);
         z = z(hi,wi);
   % z = load2P(ImageFile,'Frames',jj,'Double');
    %z = z(:,:,1,1);
  %  z=theStack(:,:,jj);

    z = z - gl(jj)*l;
    
    
    z = circshift(z,T(jj,:));
    
    z = z./thestd;
    
    
    
    z = conv2(g,g,z,'same');
    
    z = reshape(z(:) - (Q2*(Q2'*z(:))),size(z));
    
    k  =  k + (z./s2).^4;
    
end



sm = s2./m2;

k = k/length(theMatrix) - 3;


%[p,f,~] = fileparts(ImageFile{1});
%fname = fullfile(p,f);

try
%    save([fname '_align.mat'],'m','v','thestd','sm','k','T','Q','-mat');

    save(outputname,'m','v','thestd','sm','k','T','Q','-mat','-append');
    %save([fname '.alignDebug'],'msi','vsi','thestdi','gli','li','si','-mat','-append');

catch
    
    save([outputname],'m','v','thestd','sm','k','T','Q','-mat','-v7.3');
   % save([fname '.alignDebug'],'msi','vsi','thestdi','gli','li','si','-mat','-v7.3');

end


 
if computeci
fprintf(['Starting CI\n']');
sbxcomputeci3D(ImageFile,outputname,hi,wi,MD); %Takes about 10 minutes, eats up a ton of RAM
end 

%end;
