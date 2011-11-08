function [reconstruction,error] = CS2(data,percentageOfSamples)
%CS Test program for compressed sensing reconstruction - MATLAB
%OPTIMIZATION

%Initialization - data
global A;
global samples;
vec=linspace(0,1,length(data));
N=length(data);

[rs cs]=size(data);
if(rs~=1)
    data=data';
end

sound(data)

numberOfSamples=ceil(length(data)*percentageOfSamples/100);
fourier_base=[];
lim_fourier_base=[];

%Construction of measurement matrix (extracting numberOfSamples samples)
A = randn(numberOfSamples,N);
A = orth(A')';

%Construction of Fourier base
i=sqrt(-1);

K=numberOfSamples;
for j=1:N
    for k=1:N
        fourier_base(j,k) = (1/sqrt(N))*(exp((-i*2*pi/N)*(j-1)*(k-1)));
    end
    for l=1:K
        lim_fourier_base(j,l) = (1/sqrt(N))*(exp((-i*2*pi/N)*(j-1)*(k-1)));
    end
end

orig_spec=fourier_base*data';

%Application of measurement matrix and Fourier base to data -
%initialization value according to minimal energy (least squares)
samples=A*data';

init=lim_fourier_base*samples; %TODO -> check //NOT SO GOOD

A=A*fourier_base;


options = optimset('Algorithm','interior-point');
[x, fval]=fmincon(@L1Norm,init,[],[],[],[],[],[],@L2Norm,options);


reconstruction=(0.8488*pi).*(fourier_base.*exp(-1))*x; %TODO: scaling ISSUE

% subplot(3,1,1);plot(vec,real(orig_spec))
% subplot(3,1,2);plot(vec,real(x))
% subplot(3,1,3);plot(vec,reconstruction)
sound(real(reconstruction))
error=-1;
end

