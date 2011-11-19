function [reconstruction,error] = CS(data,percentageOfSamples)
%CS Test program for compressed sensing reconstruction
%TODO: Problem with non-positive definite matrices
path(path, './l1magic/Optimization');

%Initialization - data
vec=linspace(0,1,length(data));
N=length(data);

% sound(data)

numberOfSamples=ceil(length(data)*percentageOfSamples/100);
fourier_base=[numberOfSamples, numberOfSamples];

%Construction of measurement matrix (extracting numberOfSamples samples)
A = randn(numberOfSamples,N);
A = orth(A')';

%Construction of Fourier base
i=sqrt(-1);

K=numberOfSamples;
for j=1:N
    for k=1:N
        fourier_base(j,k) = (1/sqrt(N))*(exp((i*2*pi/N)*(j-1)*(k-1)));
        inv_fourier_base(j,k) = (1/sqrt(N))*(exp(-(i*2*pi/N)*(j-1)*(k-1)));
    end
end

A = A*fourier_base;

%Application of measurement matrix and Fourier base to data -
%initialization value according to minimal energy (least squares)
% A=fourier_base*A;
% Phi_i=inv_fourier_base;

samples=A*data';
init=A'*samples;

%Internal point primal-dual algorithm (Cand�s)
tic
reconstruction = l1eq_pd(init, A, [], samples, 1e-3);
toc

error=sum((abs(data'-reconstruction)).^2);

inv_fourier_base=fourier_base.*(exp(-1));

size(inv_fourier_base)
size(reconstruction)

%Plots
plot(vec,data,'b',vec,real(reconstruction),'r+')
% legend('Original data','Reconstruction using CS')
% sound(samples,round(6000/divisor))
% sound(real(reconstruction))

end
