function o = subststruct(o,s,cs)
% subststruct - Substitue tomSym symbols from a struct
%
% o = subststruct(o,s,cs) substitutes the symbols defined in struct into
% the symbolic expression o.
%
% If cs = true, then "constant" symbols are replaced by their numeric
% value.
%
% See also: subs

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

% The overloaded function actually does something. For non-tomSym input the
% output equals the input

return
