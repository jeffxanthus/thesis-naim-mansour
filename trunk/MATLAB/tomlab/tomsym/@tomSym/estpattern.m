function p = estpattern(o)
% tomSym/estpattern - Estimate the sparsity pattern of a tomSym
%
% p = estpattern(o) estimates the sparsity pattern of o, by substituting
% random values for all symbols that it contains.
%
% Although the results are usually correct, estpattern may return patterns
% that contain zeros where there should be ones. For example, the
% expression "max(0,x-0.5)" will return a random pattern containing 50 %
% zeros.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-05 by rutquist for TOMLAB release 7.7

s = symbstruct(struct,o);

sf = fieldnames(s);
for i=1:length(sf)
    s.(sf{i}) = rand(size(s.(sf{i})));
end

p = subststruct(o,s,true) ~= 0;
