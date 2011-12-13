function y = uminus(a)
% tomSym/uminus - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

%Simplify: --a = a.
if tomCmp(a,'uminus')
    y = operand(1,a);
elseif tomCmp(a,'minus')
    y = operand(2,a)-operand(1,a);
else
    y = quickop(mfilename,a);
end
