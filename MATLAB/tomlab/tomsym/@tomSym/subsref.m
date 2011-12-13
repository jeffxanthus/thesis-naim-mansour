function [o,varargout] = subsref(o,s)
% tomSym/subsref - Object index lookup

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if tomCmp(o,'cellarray')
    idx = s.subs{1};
    o = operand(idx,o);
    varargout = cell(1,nargout-1); % Why???
    return
end

if ~strcmp(s.type,'()')
    error([s.type ' indexing is not implemented tomSym']);
end

if length(s.subs)==1
    idx = s.subs{1};
    o = lookup(o,idx);
elseif length(s.subs)==2
    o = submatrix(o,s.subs{1},s.subs{2});
else
    error('A maximum of two indexes are currently supported.');
end
