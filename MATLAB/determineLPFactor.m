function [factor] = determineLPFactor(input)
%DETERMINELPFACTOR Summary of this function goes here
%   Detailed explanation goes here

[a b]=size(input);
if a~=1
    input=input';
end

N=length(input);
x=dct(input,N);

p=N/2;
if mod(N,2)~=0
    p=(N-1)/2;
end

lowerHalf=sum(abs(x(:,1:p)));
upperHalf=sum(abs(x(:,p+1:end)));

ratio=upperHalf/lowerHalf;

factor=1;
if ratio<0.02
    factor=0.25;
elseif ratio <0.05
    factor=0.55;
else
    factor=0.88;
end
end

