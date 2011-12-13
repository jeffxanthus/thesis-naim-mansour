function [x,r] = CSDeclip(data)
%CSDECLIP - data is already clipped signal
%Naim Mansour
global methodChoice
global A
global samples
global regularization


[rs cs]=size(data);
if(rs~=1)
    data=data';
end
[rs2 cs2]=size(data);
N=cs2;

M=[];
Mn=[];
Mp=[];
samples=[];
MaxS=max(data);
MinS=min(data);
%Make sure the Base decomposes the samples as part of the signal, not as
%the signal itself!
for t=1:cs2 %border values possibly added to the clipped values - SOLVED
    if(data(1,t)<MaxS && data(1,t)>MinS)
        samples=[samples data(1,t)];
    elseif(((t==1 && data(1,t+1)-data(1,t)==0)  || (t==cs2 && data(1,t)-data(1,t-1)==0))... 
            || ((t~=1 && t~=cs2) && data(1,t-1)-data(1,t+1))==0)
        M=[M t];
%         samples=[samples 0];
        if data(1,t)==MaxS
            Mp=[Mp t];
        else
            Mn=[Mn t];
        end
    else
        samples=[samples data(1,t)];
    end
end
samples=samples';

%If no clipping exists, no reason to declip
if length(M)>=1
    %Remove all rows in the DCT base, according to the sampling matrix!
    % A=DCTBase(N,N,-1);
    A=dctmtx(N)';
    B=A;
    A(M,:)=[]; 

    if length(A)==0
        disp('Clipping length exceeds frame length')
        return;
    end

    %Extra constraints
    MclA=zeros(N,N);
    MclA(Mp,:)=-B(Mp,:);
    MclA(Mn,:)=B(Mn,:);
    eps=0.9;
    offSetP=-max(samples)*eps;
    offSetN=min(samples)*eps;
    theta=(max(offSetP,offSetN)+5).*ones(N,1);
    theta(Mp,:)=offSetP;
    theta(Mn,:)=offSetN;
    

    if regularization==[]
        regularization=0.01;
    end
    if methodChoice == []
        methodChoice=3;
    end

    %Solve the constrained L1 optimization (with lambda regularization)
    switch methodChoice
        case 1
            x=SolveOMP(A,samples,N,50); %--FAST FAVORITE SO FAR
        case 2
            x=OMPDeclip(A,samples,N,MclA,theta,80); %--FAST FAVORITE SO FAR
        case 3 
            x=SolveBP(A,samples,N,50,regularization,1e-4); %Investigate parameter impact
%             options = optimset('Algorithm','interior-point','Display','on');
%             [x, fval]=fmincon(@L1Norm,offSet*ones(N,1),MclA,theta,[],[],[],[],@L2Norm,options);
        case 4
            x=IRL1(A,samples,N,50,0,1e-4); %Development in progress
    %           x=Threshold_ISD_1D(A,samples);
        case 5
    %         x=SolveLasso(A,samples,N); %--VERY SLOW, NOT THAT ACCURATE
    %           options = struct('Verbose',0);
    %           [x,niter,residuals,outputData,opts]=NESTA(A,[],samples,0.01,1e-4,options);
            disp('This option no longer exists.')
            disp('Too bad...')
            return;
    end

    r=idct(x)';


    %Magical factor - renders 2-3dB extra on the missing sample SNR
    missingRatio=length(M)/N;
    if ~(methodChoice==1 | methodChoice==2)
    if (methodChoice==1 | methodChoice==2)
        init=0.9;
       if missingRatio<=0.05
        init=0.9
       elseif missingRatio<=0.4
            init=0.8
       else
            init=0.7
       end
         limit=1;
    else
       limit=1.1;
       if missingRatio<=0.05
        limit=1.1
       elseif missingRatio<=0.4
            limit=1.2
       else
            limit=1.3
       end 
       init=1;
    end
    init
    limit
    stop=false;
    k=1;
    factor=1;
    while ~stop
        mlength=1;
        factor=1;
        while (k<=length(M)-1 && M(1,k+1)==M(1,k)+1)
            mlength=mlength+1;
            k=k+1;
        end
        mlength
        if mlength>200
            if mod(mlength,2)==0
                factor=[linspace(init,limit,mlength/2) linspace(limit,init,mlength/2)];
            else
                factor=[linspace(init,limit,(mlength-1)/2) limit*1.01 linspace(limit,init,(mlength-1)/2)];
            end
        else
            factor=ones(1,mlength);
        end
        r(1,M(1,k)-(mlength-1):M(1,k))=factor.*r(1,M(1,k)-(mlength-1):M(1,k));
        k=k+1;
        if k>length(M)
            stop=true;
        end
    end
    end
    % r=limit.*r;

    if ~(methodChoice==1 | methodChoice==2)
        data(1,M)=r(1,M);
        r=data;
    end
else
    x=dct(data);
    r=data;
end

% subplot(5,1,1);plot(data);
% title('Clipped signal')
% axis([0 N (min(data)-1) (max(data)+1)])
% subplot(5,1,2);plot(r(1:N,1));
% title('Reconstructed signal')
% axis([0 N (min(r)-1) (max(r)+1)])
% % subplot(5,1,3);plot(orig)
% % title('Original signal')
% % axis([0 N (min(orig)-1) (max(orig)+1)])
% subplot(5,1,4);plot(samples)
% title('Samples')
% axis([0 N (min(samples)-1) (max(samples)+1)])
% subplot(5,1,5);plot(abs(x))
% title('Spectral representation of reconstruction')
% axis([0 N (min(abs(x))-1) (max(abs(x))+1)])
end

