function y = isdependent(f,x)
% tomSym/isdependent - Determine if f is dependent of the symbol x

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if isa(x,'tomSym')
    x = operand(1,x);
end

if ~ischar(x)
    error('Expected a symbol name');
end

s = symbstruct(struct, f);
y = isfield(s,x);
