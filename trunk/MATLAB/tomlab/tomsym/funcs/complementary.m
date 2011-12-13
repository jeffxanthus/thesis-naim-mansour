function y = complementary(a,b)
% complementary - Test if arrays are complementary
%
% complementary(a,b) returns true if all a>=0, b>=0 and all a.*b == 0 
%
% That is: For each element i in the nonnegative arrays a and b, 
% either a(i)==0 or b(i)==0.
%
% The arrays a and b must have identical size.
%
% This function mainly exists so that it can be overloaded for tomSym,
% where it is handled efficiently by the KNITRO solver. When 
% complementary(a,b) is used as a tomSym constraint, it is mathematically
% equivalent to {a>=0, b>=0, a.*b==0}

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~all(size(a)==size(b))
    error('Complementary condition requires that arrays have the same size.');
end

y = all(a>=0) && all(b>=0) && all(a.*b==0);
