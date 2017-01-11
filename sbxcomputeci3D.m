function [] = sbxcomputeci3D(ImageFile,outputname,hi,wi,MD)

%[p,f,~] = fileparts(ImageFile{1});
%fname = fullfile(p,f);

vals = load([outputname],'-mat','T','v','Q','thestd');

T = vals.T;

Q = vals.Q;

s = sqrt(vals.v);

thestd = vals.thestd;



%Config = load2PConfig(ImageFile);



%Compute sum, sum of squares


imsize = [numel(hi) numel(wi)];
%imsize = size(vals.v); %[info.recordsPerBuffer,size(info.S,2)];

nframes = MD.hFastZ.numVolumes*length(ImageFile);
fpi = MD.hFastZ.numVolumes;


%mean, standard deviation

%compute correlation coefficient with respect to 3x3 window

winsize = 35;

res = .5;





p = gcp();

nblocks = p.NumWorkers;



xray = zeros([floor(imsize*res),winsize,winsize,nblocks],'double');

c = zeros(imsize(1),imsize(2),2,nblocks);



maxidx = nframes;

parfor ii = 1:nblocks  %nblocks % split frames evenly across cores
    
    rg = floor((ii-1)*maxidx/nblocks)+1:floor(ii*maxidx/nblocks);
    
    [c(:,:,:,ii),xray(:,:,:,:,ii)] = doOneBlock(ImageFile,imsize,res,winsize,T,Q,s,nframes,rg,thestd,MD,hi,wi);
    
end

c = sum(c,4);

xray = sum(xray,5);



c3 = (c(:,:,2)-c(:,:,1))/8/nframes;

xray = single(xray/nframes/2);

xray = int16(xray*2^15);



save(outputname,'-mat','c3','xray','-append');



end



function [c,xray] = doOneBlock(ImageFile,imsize,res,winsize,T,Q,s,nframes,rg,thestd,MD,hi,wi)

xray = zeros([floor(imsize*res),winsize,winsize],'double');

c = zeros(imsize(1),imsize(2),2);



Am = zeros([floor(imsize*res),winsize,winsize],'double');

Ar2 = 0;




for nn = rg(1):rg(length(rg))
    
    LL(:,1)=1:MD.hFastZ.numVolumes;  %index of 1: nframes per file
    LL(:,2)=1:2:MD.hFastZ.numVolumes*2; %every other frame
    LA=mod(nn,MD.hFastZ.numVolumes);  %which frame 
    if LA==0;
        LA=MD.hFastZ.numVolumes;
    end;

    A= bigread3(ImageFile{ceil(nn/MD.hFastZ.numVolumes)},LL(LA,2),1);%,1:Config(1).Frames,length(nn));

    A = double(A);
    A = A(hi,wi,:);
    A = A./thestd;                                  % normalize by ? (whiten image?)
    Ar= circshift(A,T(nn,:));                      % apply motion correction
    
    
    Ar = reshape(Ar(:) - Q*(Q'*Ar(:)),size(Ar));
    
    Ar = Ar./s;                                     % normalize by ?
    
  
    
    
    c(:,:,1) = c(:,:,1) + Ar.^2;
    
    c(:,:,2) = c(:,:,2) + conv2(ones(3,1),ones(1,3),Ar,'same').*Ar;
    
    
    
    %Ar = imresize(Ar,res);
    
    Ar = conv2([.5,1,.5],[.5,1,.5],Ar,'same')/4;    % blur
    
    Ar = Ar(2:2:end,2:2:end);                       % subsample
    
    Ar2 = Ar2 + Ar;
    
    if mod(nn,2)==0 % every other frame
        
        %At low temporal res
        
        Ar = Ar2;
        
        for ii = 1:winsize
            
            for jj = 1:winsize
                
                Am(:,:,ii,jj) = Ar.*circshift(Ar,[-(ii-ceil(winsize/2)),-(jj-ceil(winsize/2))]);
                
            end
            
        end
        
        
        
        xray = xray + Am;
        
        Ar2 = 0;
        
    end
    
    
    
    if mod(nn,100)==0
        
        fprintf('Pass 2, #%d/%d\n',nn,nframes);
        
    end
    
end

end
