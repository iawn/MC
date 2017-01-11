function [m,v,T] = sbxalignparmulti3D(ImageFile,thestd,gl,l,theMatrix,hi,wi,MD)



% Aligns images in fname for all indices in idx

%

% m - mean image after the alignment

% T - optimal translation for each frame

dep=MD.hFastZ.numFramesPerVolume;
for jj=1:dep;
    a = bigread3(ImageFile{theMatrix{jj}(1,2)},theMatrix{jj}(1,1),1);         %load2P(ImageFile,'Frames',jj,'Double');
    a = double(a);
    A{jj} = a(hi{jj},wi{jj});
    
end
p = gcp();

nblocks = 2^floor(log2(p.NumWorkers));



rg = 1:length(theMatrix);

rgs = spliteven(rg,log2(nblocks));


for jj = 1:dep
    rg1{jj} = 33:size(A{jj},1);
    
    rg2{jj} = 46:size(A{jj},2)-45;
end

for jj=1:dep;
    thestd{jj} = thestd{jj}(rg1{jj},rg2{jj});
    A{jj} = A{jj}(rg1{jj},rg2{jj});
    l{jj} = l{jj}(rg1{jj},rg2{jj});
    
    ms{jj} = zeros([size(A{jj}),nblocks]);
    vs{jj} = zeros([size(A{jj}),nblocks]);
    Ts{jj} = cell(nblocks,1);
end


for jj=1:dep;
    msi=ms{jj}; vsi=vs{jj};  Tsi=Ts{jj};
    parfor ii = 1:nblocks
        subrg = rgs{ii};
        
        
        [msi(:,:,ii),vsi(:,:,ii),Tsi{ii}] = sbxalignsub(ImageFile,subrg,rg1{jj},rg2{jj},thestd{jj},gl{jj},l{jj},theMatrix{jj},hi{jj},wi{jj});
        
        
    end;
    ms{jj}=msi; vs{jj}=vsi;  Ts{jj}=Tsi;
end;


%non-recursive
for jj=1:dep;
    for nn = 1:log2(nblocks)
        
        nblocksafter = 2^(log2(nblocks)-nn);
        
        msnew{jj} = zeros([size(A{jj}),nblocksafter]);
        
        vsnew{jj} = zeros([size(A{jj}),nblocksafter]);
        
        Tsnew{jj} = cell(nblocksafter,1);
        
        
        
        for ii = 1:nblocksafter
            
         [u{jj},v{jj}] = fftalign(ms{jj}(:,:,ii*2-1)/length(Ts{jj}{ii*2-1}), ms{jj}(:,:,ii*2  )/length(Ts{jj}{ii*2  }));
            
            
            
        Tsnew{jj}{ii} = [bsxfun(@plus,Ts{jj}{ii*2-1},[u{jj},v{jj}])
            Ts{jj}{ii*2}];
        
        Ar{jj} = circshift(ms{jj}(:,:,ii*2-1),[u{jj}, v{jj}]);
        
        msnew{jj}(:,:,ii) = (Ar{jj}+ms{jj}(:,:,ii*2  ));
        
        
        
        Ar{jj} = circshift(vs{jj}(:,:,ii*2-1),[u{jj}, v{jj}]);
        
        vsnew{jj}(:,:,ii) = (Ar{jj}+vs{jj}(:,:,ii*2  ));
        
        
        
        end
        
        
        
        
        ms{jj} = msnew{jj};
        
        vs{jj} = vsnew{jj};
        
        Ts{jj} = Tsnew{jj};
        
    end
    
    
    
    m{jj} = ms{jj}/length(theMatrix{jj});
    
    v{jj} = vs{jj}/length(theMatrix{jj});
    
    T{jj} = Ts{jj}{1};
    
end



    function A = spliteven(idx,ns)
        
        if ns > 0
            
            idx0 = idx(1:floor(end/2));
            
            idx1 = idx(floor(end/2)+1 : end);
            
            idx0s = spliteven(idx0,ns-1);
            
            idx1s = spliteven(idx1,ns-1);
            
            A = {idx0s{:},idx1s{:}};
            
        else
            
            A = {idx};
            
        end
        
        
        
    

    function [m,v,T] = sbxalignsub(ImageFile,idx,rg1,rg2,thestd,gl,l,theMatrix,hi,wi)
        
        if(length(idx)==1)
            
            %disp('L1')
            %idx
            
            
            
            A = bigread3(ImageFile{theMatrix(idx(1),2)},theMatrix(idx(1),1),1);%load2P(ImageFile,'Frames',jj,'Double');
            A = double(A);
            A = A(hi,wi);
            A = A(rg1,rg2);
            
            
            
            %       A = A(rg1,rg2);
            
            %A = A - gl(idx+1)*l;
            A = A - gl(idx)*l;
            A = A./thestd;
            
            
            
            m = A;
            
            v = A.^2;
            
            T = [0 0];
            
            
            
        elseif (length(idx)==2)
            
            
            A = bigread3(ImageFile{theMatrix(idx(1),2)},theMatrix(idx(1),1),1);%load2P(ImageFile,'Frames',jj,'Double');
            A = double(A);
            A = A(hi,wi);
            A = A(rg1,rg2);
            
            %A = A - gl(idx(1)+1)*l;
            A = A - gl(idx(1))*l;
            
            A = A./thestd;
            
            
            
            B = bigread3(ImageFile{theMatrix(idx(2),2)},theMatrix(idx(2),1),1);%load2P(ImageFile,'Frames',jj,'Double');
            B = double(B);
            B = B(hi,wi);
            B = B(rg1,rg2);
            
            %       B = B - gl(idx(2)+1)*l;
            B = B - gl(idx(2))*l;
            
            B = B./thestd;
            
            
            
            [u,v] = fftalign(A,B);
            
            
            
            Ar = circshift(A,[u,v]);
            
            m = Ar+B;
            
            T = [[u v] ; [0 0]];
            
            
            
            Ar = circshift(A.^2,[u,v]);
            
            v = (Ar+B.^2);
            
            
            
        else
            
            % disp('else')
            % idx
            
            idx0 = idx(1:floor(end/2));
            
            idx1 = idx(floor(end/2)+1 : end);
            
            [A,v1,T0] = sbxalignsub(ImageFile,idx0,rg1,rg2,thestd,gl,l,theMatrix,hi,wi);
            
            [B,v2,T1] = sbxalignsub(ImageFile,idx1,rg1,rg2,thestd,gl,l,theMatrix,hi,wi);
            
            
            
            [u,v] = fftalign(A/length(idx0), B/length(idx1));
            
            
            
            Ar = circshift(A,[u, v]);
            
            m = (Ar+B);
            
            T = [(ones(size(T0,1),1)*[u v] + T0) ; T1];
            
            
            
            v1 = circshift(v1,[u, v]);
            
            v = v1+v2;
            
            
            
        end
        
    
