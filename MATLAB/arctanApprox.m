function [ y ] = arctanApprox( x )
%UNTITLED approximation of arctan
%   Steven De Hertogh
if x >1
    y = pi/2 - 1/x + 1/(3*x^3);
else
    y = x - (x^3)/3;
end

end

