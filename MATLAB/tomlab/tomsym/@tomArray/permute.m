function o = permute(o,order)
% tomArray/permute - Overloaded function
%
%   B = permute(A,ORDER) rearranges the dimensions of A so that they are in
%   the order specified by B.
%
% If A is doesn't have named indexes, then ORDER should be a numeric
% vector, in the same manner as for Matlab's builtin permute function.
%
% If B has named indexes, then ORDER can be a cell array of index names.
%
% Example:
%
%   B = permute(A('i','j'),{'j','i'})
%
% is the same as
%
%   B = permute(A, [2 1])

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if iscell(order)
    if isempty(o.ni)
        error('Cannot use named indexes for permuting unnamed tomArray.');
    end
    if length(order)~=length(o.ni)
        error('Number of named inexes must match number of indexes in tomArray.');
    end
    if length(order)~=length(unique(order))
        error('Each index must occur exactly once.');
    end
    
    n = zeros(size(o.ni));
    for i=1:length(n)
        nn = find(strcmp(order{i}, o.ni));
        if isempty(nn)
            error(['Index ''' order{i} ''' not found in tomArray']);
        end
        n(i) = nn;
    end
else
    n = order;
end

if length(n)~=length(o.sz)
    error('ORDER must have at exactly N elements for an N-D tomArray');
end

if all(n(:)'==1:length(n))
    return;
end

idxM = reshape(1:numel(o),o.sz);
idxM = permute(idxM, n);
o.X = o.X(idxM(:));
o.sz = o.sz(n);
if ~isempty(o.ni)
    o.ni = {o.ni{n}};
end

    
