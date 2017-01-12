function [m,v,T] = sbxalignparmulti6(ImageFile,thestd,gl,l,theMatrix,hi,wi)



% Aligns images in fname for all indices in idx

% 

% m - mean image after the alignment

% T - optimal translation for each frame




 A = bigread3(ImageFile{theMatrix(1,2)},theMatrix(1,1),1);%load2P(ImageFile,'Frames',jj,'Double');
 A = double(A);
A= A(hi,wi);
 

    p = gcp();

    nblocks = 2^floor(log2(p.NumWorkers));

    

    rg = 1:length(theMatrix);
   
    rgs = spliteven(rg,log2(nblocks));    

    

    rg1 = 33:size(A,1);

    rg2 = 46:size(A,2)-45;

    

    thestd = thestd(rg1,rg2);

    

    A = A(rg1,rg2);

    l = l(rg1,rg2);

    

    ms = zeros([size(A),nblocks]);

    vs = zeros([size(A),nblocks]);

    Ts = cell(nblocks,1);

    

    parfor ii = 1:nblocks %%%parfor actively changing

        subrg = rgs{ii};

        [ms(:,:,ii),vs(:,:,ii),Ts{ii}] = sbxalignsub(ImageFile,subrg,rg1,rg2,thestd,gl,l,theMatrix,hi,wi);   

    end

    

    %non-recursive

    for nn = 1:log2(nblocks)

        nblocksafter = 2^(log2(nblocks)-nn);

        msnew = zeros([size(A),nblocksafter]);

        vsnew = zeros([size(A),nblocksafter]);

        Tsnew = cell(nblocksafter,1);

        

        for ii = 1:nblocksafter

            [u,v] = fftalign(ms(:,:,ii*2-1)/length(Ts{ii*2-1}), ms(:,:,ii*2  )/length(Ts{ii*2  }));



            Tsnew{ii} = [bsxfun(@plus,Ts{ii*2-1},[u,v]);

                         Ts{ii*2}];

            Ar = circshift(ms(:,:,ii*2-1),[u, v]);

            msnew(:,:,ii) = (Ar+ms(:,:,ii*2  ));

            

            Ar = circshift(vs(:,:,ii*2-1),[u, v]);

            vsnew(:,:,ii) = (Ar+vs(:,:,ii*2  ));



        end

        

        ms = msnew;

        vs = vsnew;

        Ts = Tsnew;

    end

    

    m = ms/length(theMatrix);

    v = vs/length(theMatrix);

    T = Ts{1};

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
%     
  
        if 0;%numel(unique(theMatrix(idx,2)))==1 %same file for each
            %read both files to save time
            dif = theMatrix(idx(2),1) - theMatrix(idx(1),1);
            [~, ~, AB] = bigread3(ImageFile{theMatrix(idx(1),2)},theMatrix(idx(1),1),dif);
            A = AB(:,:,1);
            B = AB(:,:,dif);
            
        else %if different images load seperately
            A = bigread3(ImageFile{theMatrix(idx(1),2)},theMatrix(idx(1),1),1);%load2P(ImageFile,'Frames',jj,'Double');
            B = bigread3(ImageFile{theMatrix(idx(2),2)},theMatrix(idx(2),1),1);%load2P(ImageFile,'Frames',jj,'Double');
        end
   %orig code  
%          A = bigread3(ImageFile{theMatrix(idx(1),2)},theMatrix(idx(1),1),1);%load2P(ImageFile,'Frames',jj,'Double');
         A = double(A);
         A = A(hi,wi);
         A = A(rg1,rg2);

        %A = A - gl(idx(1)+1)*l;
        A = A - gl(idx(1))*l;

        A = A./thestd;

        
%orig code
%          B = bigread3(ImageFile{theMatrix(idx(2),2)},theMatrix(idx(2),1),1);%load2P(ImageFile,'Frames',jj,'Double');
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

end
