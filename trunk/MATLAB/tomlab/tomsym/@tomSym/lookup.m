function y = lookup(v,idx)
% tomSym/lookup - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-11-21 by rutquist for TOMLAB release 7.7

if ischar(idx) && strcmp(idx,':')
    y = vec(v);
    return
end    

if islogical(idx)
    idx = find(idx);
end

if tomCmp(v,'constant')
    y = lookup(operand(1,v),idx);
    return
end

if (size(idx,1)==1 && size(idx,2)>1 && size(v,2)==1 && size(v,1)>1) || ...
    (size(idx,2)==1 && size(idx,1)>1 && size(v,1)==1 && size(v,2)>1)
    % Stupid Matlab convention.
    idx = idx';
end

if isnumeric(idx)
    % Sanity check
    if ~all(vec(idx)==round(vec(idx))) || any(vec(idx)<1)
        error('Subscript indices must either be real positive integers or logicals.');
    end
    if any(vec(idx)>numel(v))
        error('Index exceeds matrix dimensions.');
    end
end

% Eliminate repeated lookup
if tomCmp(v,'lookup')
    y = lookup(operand(1,v),lookup(operand(2,v),idx));
    return
end

% Eliminate setSymmetric
if tomCmp(v,'setSymmetric') && isnumeric(idx)
    n = size(v,1);
    it = setSymmetric(1:n*(n+1)/2);
    y = lookup(operand(1,v),it(idx));
    return
end

% Eliminate reshape
if tomCmp(v,'reshape') && isnumeric(idx)
    y = lookup(operand(1,v),idx);
    return
end

% Possibly eliminate transpose
if tomCmp(v,'ctranspose') && size(v,1)~=1 && size(v,2)~=1
    vn = reshape(1:numel(v),size(v,2),size(v,1));
    it = lookup(vn',idx);
    y = lookup(operand(1,v),it);
    return
end
if tomCmp(v,'ctranspose') && (size(v,1)==1 || size(v,2)==1) && numel(idx)==1;
    y = lookup(operand(1,v),idx);
    return
end

% Possibly eliminate setdiag
if tomCmp(v,'setdiag') && isnumeric(idx)
    [idxi, idxj] = ind2sub(size(v),idx);
    if all(idxi==idxj)
        y = lookup(operand(1,v),idxi);
        return
    end
end

% Possibly eliminate sparse
if tomCmp(v,'sparse') && isnumeric(idx) && isnumeric(operand(1,v)) && isnumeric(operand(2,v))
    idxm = sparse(operand(1,v),operand(2,v),1:length(operand(1,v)),size(v,1),size(v,2));
    idxs = lookup(idxm,idx);
    idxf = find(idxs);
    if isempty(idxf)
        y = spzeros(size(idx));
    elseif length(idxf)==length(idxs)
        y = reshape(lookup(operand(3,v),idxs),size(idx));
    else
        [ii,jj] = ind2sub(size(idx),idxf);
        y = sparse(ii,jj,lookup(operand(3,v),idxs(idxf)),size(idx,1),size(idx,2));
    end
    return
end

% Eliminate repmat
if tomCmp(v,'repmat')
    vv = operand(1,v);
    if numel(vv)==1
        % All elements are identical. Lookup just becomes another repmat.
        y = repmat(vv,size(idx,1),size(idx,2));
    else
        iv = reshape(1:numel(vv),size(vv));
        ii = repmat(iv,operand(2,v),operand(3,v));
        y = lookup(vv,lookup(ii,idx));
        if size(vv,2)==1 && size(v,2)>1 && size(idx,1)==1;
            % Stupid Matlab transposing vectors arbitrarily
            y = y';
        end
    end
    return
end

if numel(idx) < numel(v)
    % Simplifications that make sense when the result is smaller than
    % the original v.

    if isa(v,'tomSym') && any(strcmp(operator(v), {'uminus','uplus', ...
            'cos','acos','log','exp','cot','sec','csc',...
            'acot','asec','acsc','cosh','acosh','acoth','asech','acsch',...
            'log2','log10','sin','asin','tan','atan','sinh','asinh',...
            'tanh','atanh','abs','sign'}))
        % Unary element-wise functions
        y = feval(operator(v), lookup(operand(1,v),idx));
        return
    elseif isa(v,'tomSym') && any(strcmp(operator(v), {'plus','minus', ...
            'times','rdivide','gt','lt','ne','ge','le','eq'}))
        % Binary element-wise functions
        y = feval(operator(v), lookup(operand(1,v),idx), lookup(operand(2,v),idx));
        return
    elseif tomCmp(v,'smplus')
        y = operand(1,v)+lookup(operand(2,v),idx);
        return
    elseif tomCmp(v,'smtimes')
        y = operand(1,v)*lookup(operand(2,v),idx);
        return
    elseif tomCmp(v,'mtimes') && numel(idx)==1 && isnumeric(idx)
        [i,j] = ind2sub(size(v),idx);
        y = submatrix(operand(1,v),i,':')*submatrix(operand(2,v),':',j);
        return
    elseif tomCmp(v,'spower')
        y = lookup(operand(1,v),idx).^operand(2,v);
        return
    elseif tomCmp(v,'vecsubs') && size(v,1)==1 && numel(idx)==1
        vo = operands(v);
        vo{2} = lookup(vo{2},idx);
        y = vecsubs(vo{:});
        return;
    elseif tomCmp(v,'vecsubs') && numel(idx)==1 && isnumeric(idx)
        vo = operands(v);
        sz = size(vo{1}.fun);
        p = ceil(idx/(sz(1)*sz(2)));
        vo{2} = lookup(vo{2},p);
        y = lookup(vecsubs(vo{:}),idx-(p-1)*sz(1)*sz(2));
        return;
    elseif tomCmp(v,'submatrix')
        if size(v,1)==1
            y = submatrix(v,':',idx);
            return
        elseif size(v,2)==1
            y = submatrix(v,idx,':');
            return
        end
    end
end
    
% Possibly eliminate concatenations
if isnumeric(idx) && (tomCmp(v,'horzcat') || tomCmp(v,'vertcat'))
    isvert = tomCmp(v,'vertcat');
    iv = zeros(size(v));
    ix = 0;
    if isvert
        for i = 1:nOperands(v)
            iv(ix+1:ix+size(operand(i,v),1),:) = i;
            ix = ix+size(operand(i,v),1);
        end
    else
        for i = 1:nOperands(v)
            iv(:,ix+1:ix+size(operand(i,v),2)) = i;
            ix = ix+size(operand(i,v),2);
        end
    end
    lv = lookup(iv,idx);
    us = false(1,nOperands(v));
    for i=1:length(us)
        us(i) = any(lv(:)==i);
    end
    if ~all(us)
        % Do elimination
        fu = find(us);
        cv = cell(1,length(fu));
        iv = false(size(v));
        ix = 0;
        if isvert
            for i = 1:nOperands(v)
                iv(ix+1:ix+size(operand(i,v),1),:) = us(i);
                ix = ix+size(operand(i,v),1);
            end
        else
            for i = 1:nOperands(v)
                iv(:,ix+1:ix+size(operand(i,v),2)) = us(i);
                ix = ix+size(operand(i,v),2);
            end
        end
        iv = find(iv);
        ri = zeros(size(v));
        ri(iv) = 1:length(iv);
        for i = 1:length(fu)
            cv{i} = operand(fu(i),v);
        end
        if isvert
            v = vertcat(cv{:});
        else
            v = horzcat(cv{:});
        end
        idx = lookup(ri,idx);
        y = lookup(v,idx);
        return
    end
end

if numel(idx)==numel(v) && isnumeric(idx) && all(idx(:)==(1:numel(idx))')
    % Index covers all elements of v
    y = reshape(v,size(idx));
elseif numel(v) == 1
    % Lookup into a single-element list
    y = repmat(v,size(idx));
elseif isnumeric(v) && (numel(v)==2 || v(2)==0.5*(v(1)+v(3)) && nnz(diff(diff(v(:))))==0)
    % V is a a linear mapping
    k = v(2)-v(1);
    m = v(1)-k;
    y = k*idx+m;
else
    y = tomSym(mfilename,size(idx,1),size(idx,2),v,idx);
end

