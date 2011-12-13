function y = constpart(o)
% tomSym/constpart - The zeroth term of the McLaurin series of a tomSym
%
% y = constpart(o) evaluates the object o with all symbols set equal to
% zero, leaving only the constant term.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-05 by rutquist for TOMLAB release 7.7

symbs = symbols(o,'struct');

s = struct;
sl = fieldnames(symbs);
for i=1:length(sl)
    s.(sl{i}) = zeros(size(symbs.(sl{i})));
end

y = subststruct(o,s,true);
