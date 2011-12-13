function y = positiveSemidefinite(a)
% tomSym/positiveSemidefinite - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if size(a,1)~=size(a,2)
    error('Matrix must be square.');
end

% Simplify for scalar case
if size(a,1)==1
    if tomCmp(a,'uminus')
        y = -a<=0;
    else
        y = a>=0;
    end
else
    y = tomSym(mfilename,1,1,a);
end
