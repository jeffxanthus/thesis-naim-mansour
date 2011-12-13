function o = repmat(o,m,n)
% tomArray/repmat - Overloaded function.
%
% B = repmat(A,M,N) creates M-by-N tiling copies of A.
%
% B = repmat(A,[M N]) is the same as repmat(A,M,N).
% 
% B = repmat(A,[M N P ...]) works in any number of dimensions.
%
% B = repmat(A,'k',N) for a tomArray with named indexes, duplicates the
% array along the dimension 'k'. If 'k' is not one of the named indexes of
% A, then it is added to the end of the index list, increasing the number
% of dimensions by one.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    error('Requires at least 2 inputs.');
elseif nargin == 2
    siz = m;
else % nargin == 3
    if isnumeric(m)
        siz = [m n];
    elseif ischar(m)
        % Named index
        if ~isvarname(m)
            error(['Illegal index: ' s.subs{i}]);
        end

        if isempty(o.ni)
            error('Repmat with named index only works if tomArray indexes are already named.');
        end
        siz = ones(size(o.sz));
        ix = find(strcmp(o.ni,m));
        if isempty(ix)
            siz(end+1) = n;
            o.ni{end+1} = m;
        else
            siz(ix) = n;
        end
    else
        error(['Illegal index type: ' class(m)]);
    end
end

if length(o.sz)<length(siz)
    o.sz(end+1:length(siz)) = 1;
end

lst = find(o.sz>1); lst = lst(end); %find(...,'last') is incompatible with matlab6.5
fst = find(siz>1); fst = fst(1);

if isempty(fst)
    return;
end

if isempty(lst) || lst<=fst
    % Duplication only along the last dimension(s)
    o.X = repmat(o.X(:),prod(siz),1);
else
    % Duplication along earlier dimension
    idxM = reshape(1:numel(o),o.sz);
    idxM = repmat(idxM,siz);
    o.X = o.X(idxM);
end

o.sz = o.sz.*siz;

checkIndexes(o);
