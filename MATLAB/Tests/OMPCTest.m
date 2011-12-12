function [ output_args ] = OMPCTest( input_args )
%OMPCTEST 

signal=sin(linspace(0,4*pi,100));
B=dctmtx(length(signal))';
xx=dct(signal);

data=Clip(signal,0.6);

M=[];
Mn=[];
Mp=[];
samples=[];
MaxS=max(data);
MinS=min(data);

[rs2 cs2]=size(data);
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

N=cs2
%Extra constraints
    MclA=zeros(N,N);
    MclA(Mp,:)=-B(Mp,:);
    MclA(Mn,:)=B(Mn,:);
    eps=0.9;
    offSet=-max(abs(samples))*eps;
    theta=Inf.*ones(N,1);
    theta(Mp,:)=offSet;
    theta(Mn,:)=offSet;
    Mp 
    Mn
    pause
    MclA
    theta
    plot(data)
    pause
    size(xx)
    size(B)
    plot(MclA*xx','.')
    size(MclA*xx')
    if MclA*xx'-theta<=zeros(N,1)
        disp('Success')
    end

end

