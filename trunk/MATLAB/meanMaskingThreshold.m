function [ threshold ] = meanMaskingThreshold( input )
%meanMaskingThreshold Calculates mean maskingthreshold over whole input
%Steven De Hertogh

N = 512;

numberofframes = floor(length(input)/N);
masking = maskingThreshold(input);
for n = 2:numberofframes
   masking = (masking + maskingThreshold(input(N*(n-1):N*(n))))./2;
end

threshold = masking;
figure();
plot(threshold);
end

