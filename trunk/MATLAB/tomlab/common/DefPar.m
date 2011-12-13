% function x=DefPar(s,f,xdef)
%
% If 'f' is a field in structure 's' and s.f is not empty, return it.
% If 'f' is not a field in s, or s.f is empty, return xdef.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.7.0 $
% Written Dec 1, 2003.    Last modified Mar 22, 2005.

function x=DefPar(s,f,xdef)

if(nargin < 3)
    xdef = [];
end

if isfield(s,f)
    x = getfield(s,f);
    if isempty(x), x = xdef; end
else
    x = xdef;
end