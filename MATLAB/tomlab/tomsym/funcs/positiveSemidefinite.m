function y = positiveSemidefinite(A,margin)
% positiveSemidefinite - Test if a matrix is positive definite
%
% positiveSemidefinite(A) returns true if A is positive semidefinite.
%
% A must be real-valued and square.
% If it is not symmetric then 1/2*(A+A') is used.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    margin = 0;
end

if ~isreal(A)
    error('Not implemented for complex numbers yet.')
end
if size(A,1)~=size(A,2)
    error('Matrix must be square.');
end

A = 0.5*(A+A');

opts = struct('issym',true,'isreal',true,'disp',0);

if size(A,1)>1
    y = eigs(A,1,'SA',opts)>-margin;
else
    y = A>=0;
end
