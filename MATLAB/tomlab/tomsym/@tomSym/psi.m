function y = psi(k,a)
% tomSym/psi - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin==1
    % Default k=0
    a = k;
    k = 0;
end

if tomCmp(a,'ctranspose')
    % Move out transpose operator
    y = tomSym(mfilename,size(a,2),size(a,1),k,operand(1,a))';
else
    y = tomSym(mfilename,size(a,1),size(a,2),k,a);
end
