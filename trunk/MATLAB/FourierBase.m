function [F] = FourierBase(r,c,inv)
%FOURIERBASE
% N=c, since r <= c (underdetermined)
%inv: parameter that allows choice between normal or inverse Fourier transform

if(~(inv==1||inv==-1))
    disp('Incorrect argument')
    return
end

N=c;

F=zeros(r,c);

for k=1:r
F(k,:)=exp(-2*pi*i*inv*(k-1)*(0:(c-1))/N); 
end

F=F./sqrt(N);
end
