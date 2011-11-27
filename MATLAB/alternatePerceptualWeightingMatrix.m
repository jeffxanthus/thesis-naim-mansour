function [ H ] = alternatePerceptualWeightingMatrix( maskingThreshold )
%UNTITLED Summary of this function goes here
%   Steven De Hertogh

K = length(maskingThreshold);
A = 1./maskingThreshold;

H = [];
for i = 1:K
    H(i,i) = abs(A(i));
end


end

