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



parfor ii = 1:max_idx

%     A = load2P(ImageFile,'Frames',ii);
%   %  A = theStack(:,:,ii);
%     A=A(:,:,1,1);
%     
    
         A = bigread3(ImageFile{theMatrix(ii,2)},theMatrix(ii,1),1);%load2P(ImageFile,'Frames',jj,'Double');
         A = double(A);
    
         A = A(hi,wi);
    
    
    A = double(A(rg1,rg2));

    A = A - gl(ii)*l;

    A = A./thestd;

    [dx,dy] = fftalign(A,m0);

    T(ii,:) = [dx,dy];

    

    Ar = circshift(A,[dx, dy]);

    m = m+double(Ar);

    

    Ar = circshift(A.^2,[dx, dy]);

    v = v + double(Ar);



end



m = m/max_idx;

v = m/max_idx;
