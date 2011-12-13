function y=plus(a,b)
% tomSym/plus - Overload the plus operator

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-08 by rutquist for TOMLAB release 7.7

if numel(a)==1 && numel(b)~=1
    y = smplus(a,b);
elseif numel(b)==1 && numel(a)~=1
    y = smplus(b,a);
elseif tomCmp(a,'repmat') && numel(operand(1,a))==1
    y = operand(1,a)+b;
elseif tomCmp(b,'repmat') && numel(operand(1,b))==1
    y = a+operand(1,b);
elseif iszero(a)
    y = b;
elseif iszero(b)
    y = a;
elseif isnumeric(a) && numel(a)>1 && all(size(a)==size(b)) && nnz(diff(a(:)))==0
    y = smplus(a(1),b);
elseif isnumeric(b) && numel(b)>1 && all(size(a)==size(b)) && nnz(diff(b(:)))==0
    y = smplus(b(1),a);
elseif tomCmp(a,'setrows') && tomCmp(b,'setrows')
    oa = operands(a);
    ob = operands(b);
    y = setrows(oa{:},ob{2:end});
elseif (tomCmp(a,'sparse') || tomCmp(a,'setdiag')) && ...
        (tomCmp(b,'sparse') || tomCmp(b,'setdiag'))
    % plus can be used to combine sparse matrices
    if ~all(size(a)==size(b))
        error('Matrix dimensions must agree');
    end
    if tomCmp(a,'setdiag')
        ia = (1:size(a,1))';
        ja = (1:size(a,1))';
        ca = vec(operand(1,a));
    else
        ia = operand(1,a);
        ja = operand(2,a);
        ca = operand(3,a);
    end
    if tomCmp(b,'setdiag')
        ib = (1:size(b,1))';
        jb = (1:size(b,1))';
        cb = vec(operand(1,b));
    else
        ib = operand(1,b);
        jb = operand(2,b);
        cb = operand(3,b);
    end
    if isequal(ia,ib) && isequal(ja,jb)
        % Move plus inside sparse
        y = sparse(ia,ja,ca+cb,size(a,1),size(a,2));
    else
        if isnumeric(ia) && isnumeric(ja) && isnumeric(ib) && isnumeric(jb) && ...
                isempty(intersect([ia ja],[ib, jb],'rows'));
            % Combine sparse matrices
            y = sparse([ia; ib],[ja; jb],...
                [vec(ca);vec(cb)],size(a,1),size(a,2));
        else
            % Nothing to be done
            y = quickop(mfilename,a,b);            
        end
    end
elseif tomCmp(b,'uminus')
    y = a-operand(1,b);
elseif tomCmp(a,'uminus')
    y = b-operand(1,a);
elseif tomCmp(a,'minus') && isequal(operand(2,a),b)
    y = operand(1,a);
else
    y = quickop(mfilename,a,b);
end
