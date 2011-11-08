function [D] = DCTBase(r,c,inv)
%DCTBASE - implements DCT-II base of dimensions rxc
%Naim Mansour
if(~(inv==1||inv==-1))
    disp('Incorrect argument')
    return
end

N=c;
D=zeros(r,c);

if inv==1
    D(1,:)=(1/sqrt(N)).*cos((pi/N)*(0:(c-1)+0.5));
    for k=2:r
    D(k,:)=cos((pi/N)*(0:(c-1)+0.5)*k);
    end
end
if inv==-1
    D(1,:)=(1/sqrt(N)).*cos((pi/N)*(0:(c-1)+0.5));
    for k=2:r
    D(k,:)=cos((pi/N)*(0:(c-1)+0.5)*k);
    end
end

D=D.*sqrt(2/N);
end

