function y = logm(a)
% tomSym/logm - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

% Simplifications
if numel(a)==1
    y = log(a);
elseif isallone(a)
    y = zeros(size(a));
elseif strcmp(operator(a),'expm')
    y = operand(1,a);
else
    y = quickop(mfilename,a);
end
