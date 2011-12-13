function o = subsref(o,s)
% tomArray/subsref - Lookup and/or name indexes.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-26 by rutquist for TOMLAB release 7.7

if ~strcmp(s.type,'()')
    error([s.type ' indexing is not implemented for tomArray']);
end

if length(s.subs)==1 && length(o.sz)>1 && ~ischar(s.subs{1})
    idx = s.subs{1};
    sz = size(idx);
    if prod(sz)==max(sz)
        sz = max(sz);
    end
    o = tomArray(lookup(o.X,idx(:)),sz);
elseif length(s.subs)==1 && ischar(s.subs{1}) && strcmp(s.subs{1},':')
    o = o.X(:);
elseif length(s.subs)==length(o.sz)
    % Submatrix and/or naming operation
    idx = cell(size(o.sz));
    isfull = true;
    for i=1:length(o.sz)
        if ischar(s.subs{i})
            c = s.subs{i};
            idx{i} = 1:o.sz(i);
        elseif iscell(s.subs{i})
            c = s.subs{i}{1};
            idx{i} = s.subs{i}{2}(:)';
        elseif isa(s.subs{i},'tomArrayIdx')
            c = char(s.subs{i});
            idx{i} = double(s.subs{i});
            if ischar(idx{i})
                idx{i} = 1:o.sz(i);
            end
        else
            c = '';
            idx{i} = s.subs{i}(:)';
        end
        if ~( numel(idx{i}) == o.sz(i) && all(idx{i}==1:o.sz(i)) )
            isfull = false;
        end
        if islogical(idx)
            idx = find(idx);
        end
        if any(idx{i}<1) || any(idx{i}~=round(idx{i}))
            error(['TomArray index ' c ' is not a positive integers or a logical.']);
        end
        if any(idx{i}>o.sz(i))
            error(['TomArray index ' c ' exceeds matrix dimensions.']);
        end
        if ~isempty(c)
            if isvarname(c)
                if length(o.ni) >= i && ~strcmp(o.ni{i}, c)
                    error('Cannot rename index. (Tip: Use permute() to change indexing order.)');
                end
                o.ni{i} = c;
            elseif ~strcmp(c,':')
                error(['Illegal index: ' s.subs{i}]);
            end
        end
    end
    
    % Check for diagonal extraction
    rix = 1:length(o.sz);
    for i=2:length(o.ni)
        for j=1:i-1
            if ~isempty(o.ni{i}) && strcmp(o.ni{i},o.ni{j});
                rix(i) = j;
            end
        end
    end
    keep = rix==1:length(rix);
    if ~all(keep)
        % extract diagonals
        for i=1:length(rix)
            if length(idx{i}) ~= length(idx{rix(i)})
                error('Diagonal index length missmatch.');
            end
        end
        kix = find(keep);
        ix0 = cell(1,length(kix));
        for i=1:length(ix0)
            ix0{i} = 1:length(idx{kix(i)});
        end
        ix1 = cell(1,length(kix));
        if length(ix0)>1
            [ix1{:}] = ndgrid(ix0{:});
        else
            ix1 = ix0; % Stupid ndgrid can't handle single argument correctly.
        end
        ix = cell(1,length(o.sz));
        rrix = cumsum(keep);
        for i=1:length(ix)
            ix{i} = vec(idx{i}(ix1{rrix(rix(i))}));
        end        
        o.X = o.X(sub2ind(o.sz,ix{:}));
        o.sz = zeros(1,length(kix));
        for i=1:length(o.sz)
            o.sz(i) = length(idx{kix(i)});
        end
        o.ni = o.ni(keep(1:length(o.ni)));
    elseif ~isfull
        % Submatrix extraction
        idxM = getIdxM(o.sz);
        idxM = idxM(idx{:});
        if length(idx) > 1
            o.sz = size(idxM);
        else
            o.sz = length(idx{1});
        end
        o.X = o.X(idxM(:));
    end
    
    if ~isempty(o.ni);
        o = squeeze(o);
        checkIndexes(o);
    end
else
    error('Wrong number of indexes for tomArray.');
end
