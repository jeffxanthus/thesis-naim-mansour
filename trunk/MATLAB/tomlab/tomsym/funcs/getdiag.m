function v = getdiag(M)
% getdiag - Get the main diagonal from a square matrix
%
% v = getdiag(M) extracts the main diagonal from the matrix M.
%
% While getdiag does approximately the same thing as diag, it is more clear
% what it does, since diag also does the job of setdiag.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if size(M,1)~=size(M,2)
    error('Getdiag is only defined for square matrices');
end

v = diag(M);
