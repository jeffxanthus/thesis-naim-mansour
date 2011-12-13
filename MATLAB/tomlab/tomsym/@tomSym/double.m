function y = double(o)
% tomSym/double - Convert tomSym constant to double

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-08-22 by rutquist for TOMLAB release 7.7

y = subststruct(o,[],1);
    
if ~isnumeric(y);
    error('Cannot convert symbolic tomSym object to double.');
end

y = double(y);
