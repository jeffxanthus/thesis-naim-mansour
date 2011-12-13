function [ H ] = alternatePerceptualWeightingMatrix( maskingThreshold )
%UNTITLED Summary of this function goes here
%   Steven De Hertogh

K = length(maskingThreshold);
A = 10.^((1./maskingThreshold)/20);

H = [];
for i = 1:K
    H(i,i) = A(i);
    if i > K-8
        H(i,i) = A(K-8);
    end
end


end

