function y = submatrix(M,i,j)
% tomSym/submatrix - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-26 by rutquist for TOMLAB release 7.7

%if tomCmp(M,'constant')
%    y = submatrix(operand(1,M),i,j);
%    return
%end

if ischar(i) && strcmp(i,':')
    ni = size(M,1);
    ai = true;
else
    if islogical(i)
        i = find(i);
    end
    ni = length(i);
    ai = (ni==size(M,1) && isnumeric(i) && all(i(:)==(1:size(M,1))'));
end

if ischar(j) && strcmp(j,':')
    nj = size(M,2);
    aj = true;
else
    if islogical(j)
        j = find(j);
    end
    nj = length(j);
    aj = (nj==size(M,2) && isnumeric(j) && all(j(:)==(1:size(M,2))'));
end

if (isnumeric(i) && (any(i<1) || any(i>size(M,1)))) || ...
        (isnumeric(j) && ( any(j<1) || any(j>size(M,2))))
    error('Index out of range');
end

% Check if empty index. (Return empty matrix.)
if isempty(i)
    y = spzeros(0,nj);
    return
elseif isempty(j)
    y = spzeros(ni,0);
    return
end

if ~(ai || size(i,1)==1 || size(i,2)==1) || ...
        ~(aj || size(j,1)==1 || size(j,2)==1)
    error('Index lists must be vectors');
end

if size(M,1)==2 && isa(i,'tomSym')
    % Two points - use linear interpolation
    % (This makes it possible to compute derivatives for Linear
    % Programming.)
    M1 = submatrix(M,1,j);
    M2 = submatrix(M,2,j);
    p = vec(i)-1;
    y = (1-p)*M1 + p*M2;
    return
elseif size(M,2)==2 && isa(j,'tomSym')
    M1 = submatrix(M,i,1);
    M2 = submatrix(M,i,2);
    p = vec(j)'-1;
    y = M1*(1-p) + M2*p;
    return
end

if tomCmp(M,'submatrix')
    % Submatrix of submatrix
    ii = operand(2,M);
    if ~(ischar(ii) && strcmp(ii,':'))
        i = lookup(ii,i);
    end
    jj = operand(3,M);
    if ~(ischar(jj) && strcmp(jj,':'))
        j = lookup(jj,j);
    end
    y = submatrix(operand(1,M),i,j);
    return
end
    
% Do lookup instead, if just one element
if isnumeric(i) && ni==1 && isnumeric(j) && nj==1
    y = lookup(M,sub2ind(size(M),i,j));
    return
end

% Return whole M if M(:,:)
% Do lookup instead, if vector
if ai && aj
    y = M;
    return
elseif size(M,1)==1 && ai
    y = reshape(lookup(M,j),ni,nj);
    return
elseif size(M,2)==1 && aj
    y = reshape(lookup(M,i),ni,nj);
    return
end

if ni*nj < numel(M)
    % Simplifications that make sense when the result is smaller than
    % the original M.

    if isa(M,'tomSym') && any(strcmp(operator(M), {'uminus','uplus', ...
            'cos','acos','log','exp','cot','sec','csc',...
            'acot','asec','acsc','cosh','acosh','acoth','asech','acsch',...
            'log2','log10','sin','asin','tan','atan','sinh','asinh',...
            'tanh','atanh','abs','sign'}))
        % Unary element-wise functions
        y = feval(operator(M), submatrix(operand(1,M),i,j));
        return
    elseif isa(M,'tomSym') && any(strcmp(operator(M), {'plus','minus', ...
            'times','rdivide','gt','lt','ne','ge','le','eq'}))
        % Binary element-wise functions
        y = builtin('feval', operator(M), submatrix(operand(1,M),i,j), submatrix(operand(2,M),i,j));
        return
    elseif tomCmp(M,'smplus')
        y = operand(1,M)+submatrix(operand(2,M),i,j);
        return
    elseif tomCmp(M,'smtimes')
        y = operand(1,M)*submatrix(operand(2,M),i,j);
        return
    elseif tomCmp(M,'mtimes')
        y = submatrix(operand(1,M),i,':')*submatrix(operand(2,M),':',j);
        return
    elseif tomCmp(M,'sparse')
        ii = operand(1,M);
        jj = operand(2,M);
        if isnumeric(i) && isnumeric(j)
            iic = zeros(1,size(M,1));
            jjc = zeros(1,size(M,2));
            iic(i) = 1:length(i);
            jjc(j) = 1:length(j);
            u = find(iic(ii) & jjc(jj));
            y = sparse(lookup(iic,lookup(ii,u)), lookup(jjc,lookup(jj,u)), ...
                lookup(operand(3,M),u), ni, nj);
            return
        elseif ai && isnumeric(j)
            jjc = zeros(1,size(M,2));
            jjc(j) = 1:length(j);
            u = find(jjc(jj));
            y = sparse(lookup(ii,u), lookup(jjc,lookup(jj,u)), ...
                lookup(operand(3,M),u), ni, nj);
            return
        elseif isnumeric(i) && aj
            iic = zeros(1,size(M,1));
            iic(i) = 1:length(i);
            u = find(iic(ii));
            y = sparse(lookup(iic,lookup(ii,u)), lookup(jj,u), ...
                lookup(operand(3,M),u), ni, nj);
            return
        end
    elseif tomCmp(M,'vecsubs')        
        vo = operands(M);
        if size(vo{1}.fun,2)==1
            % TODO: Simplify for size(vo{1}.fun,2)>1
            if ~ai
                vo{1}.fun = lookup(vo{1}.fun,i);
            end
            vo{2} = lookup(vo{2},j);
            y = vecsubs(vo{:});
            return
        end
    elseif tomCmp(M,'kron')
        y = subkron(operand(1,M),operand(2,M),i,j);
        return
    end
end

% Possibly eliminate concatenations
if (isnumeric(i) && tomCmp(M,'vertcat')) || (isnumeric(j) && tomCmp(M,'horzcat'))
    if tomCmp(M,'vertcat');
        ij = i;
        szi = 1;
    else
        ij = j;
        szi = 2;
    end
    iv = zeros(size(M,szi));
    ix = 0;
    for k = 1:nOperands(M)
        iv(ix+1:ix+size(operand(k,M),szi)) = k;
        ix = ix+size(operand(k,M),szi);
    end
    lv = lookup(iv,ij);
    us = false(1,nOperands(M));
    us(lv) = true;    
    if ~all(us)
        % Do elimination
        fu = find(us);
        cv = cell(1,length(fu));
        iv = false(size(M,szi),1);
        ix = 0;
        for k = 1:nOperands(M)
            iv(ix+1:ix+size(operand(k,M),szi)) = us(k);
            ix = ix+size(operand(k,M),szi);
        end
        iv = find(iv);
        ri = zeros(size(M,szi),1);
        ri(iv) = 1:length(iv);
        for k = 1:length(fu)
            cv{k} = operand(fu(k),M);
        end
        ij = lookup(ri,ij);
        if szi==1
            M = vertcat(cv{:});
            y = submatrix(M,ij(:)',j);
        else
            M = horzcat(cv{:});
            y = submatrix(M,i,ij(:)');
        end
        return
    end
end

if tomCmp(M,'ctranspose')
    % Move out transpose operator
    y = submatrix(operand(1,M),j,i)';
else
    y = tomSym(mfilename,ni,nj,M,i,j);
end
