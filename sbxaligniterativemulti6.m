function [m,v,T] = sbxaligniterativemulti6(ImageFile,m0,rg1,rg2,thestd,gl,l,theMatrix,hi,wi)

% Aligns images in fname for all indices in idx

% 

% m - mean image after the alignment

% T - optimal translation for each frame




%Config = load2PConfig(ImageFile);
%fpc=Config.Frames;
%lConfig=length(Config);
%Config=Config(1);
%Config.Frames = fpc*lConfig;

T = zeros(length(theMatrix),2);

m = zeros(length(rg1),length(rg2));

v = zeros(length(rg1),length(rg2));



l = l(rg1,rg2);



max_idx = length(theMatrix);

nFile = numel(ImageFile);

if nFile>1
    nFrames = find(theMatrix(:,2)>1,1)-1;
else
    nFrames = size(theMatrix,1);
end

if nFrames>500
    fprintf('Warning large individual file sizes. maybe you should break this into chunks\n')
end

%nChunk=1; %if you add code to chunk this will change
range = theMatrix(nFrames,1)-theMatrix(1,1)+1;



for i=1:nFile

    [~,~,A1] = bigread3(ImageFile{i},theMatrix(1,1),range); 
    A1 = A1(:,:,theMatrix(1:nFrames,1)-theMatrix(1,1)+1);
    
    glReshape = reshape(gl,[numel(gl)/nFile nFile]); %so you don't have to send all of gl
parfor ii = 1:nFrames 

%     A = load2P(ImageFile,'Frames',ii);
%   %  A = theStack(:,:,ii);
%     A=A(:,:,1,1);
%     
    
         %A = bigread3(ImageFile{theMatrix(ii,2)},theMatrix(ii,1),1);%load2P(ImageFile,'Frames',jj,'Double');
         
         A = double(A1(:,:,ii));
    
         A = A(hi,wi);
    
    
    A = double(A(rg1,rg2));

    A = A - glReshape(ii,i)*l;

    A = A./thestd;

    [dx,dy] = fftalign(A,m0);

    Tpart(ii,:) = [dx,dy]; 

    

    Ar = circshift(A,[dx, dy]);

    m = m+double(Ar);

    

    Ar = circshift(A.^2,[dx, dy]);

    v = v + double(Ar);



end
T((i-1)*nFrames+1:i*nFrames,:) = Tpart; %process each chunk of the translation matrix and then concatenate together

end


m = m/max_idx;

v = m/max_idx;
