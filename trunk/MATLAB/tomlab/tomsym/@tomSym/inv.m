function y = inv(a)
% tomSym/inv - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if size(a,1)~=size(a,2)
    error('Matrix must be square');
end

if isone(a)
    % Inverse of identity matrix
    y = a;
elseif tomCmp(a,'inv')
    % Inverse of an inverse
    y = operand(1,a);
elseif tomCmp(a,'setdiag')
    % Diagonal matrix is inverted element-wise
    y = setdiag(1./operand(1,a));
else
    y = quickop(mfilename,a);
end
