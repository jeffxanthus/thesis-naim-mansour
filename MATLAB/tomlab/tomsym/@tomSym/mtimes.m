function y=mtimes(a,b)
% tomSym/mtimes - Overload the multiplication operator
%
% The * operator can mean either matrix-matrix, scalar-matrix, or
% matrix-scalar multiplication. We map them as follows:
%   mtimes  - matrix-matrix multiplication, including the case where a is
%             m-by-1 and b is 1-by-1, or a is 1-by-1 and b is 1-by-n.
%             (In particular, this includes scalar-scalar multiplication.)
%   smtimes - scalar-matrix multiplication, where a is a scalar. It is
%             assumed that this operator commutes, so matrix-scalar
%             multiplication is mapped to this, with the order of operands
%             changed.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-09-08 by rutquist for TOMLAB release 7.7

if islogical(a)
    a = double(a);
end
if islogical(b)
    b = double(b);
end

if ~(isnumeric(a) || isa(a,'tomSym'))
   error(['Argument 1 is a ' class(a) '. Arguments for tomSym/mtimes must be either symbolic or numeric.']); 
end
if ~(isnumeric(b) || isa(b,'tomSym'))
   error(['Argument 2 is a ' class(b) '. Arguments for tomSym/mtimes must be either symbolic or numeric.']); 
end

% Check for complex constants
if isnumeric(a) && ~isreal(a) && any(imag(a(:)))
    if ~isequal(a,1j);
        warning('tomsym:complexdata','Complex data encountered. Attempting to convert to tomCmplx.');
    end
    b = tomCmplx(b,spzeros(size(b)));
    y = feval('mtimes',a,b);
    return
end
if isnumeric(b) && ~isreal(b) && any(imag(b(:)))
    if ~isequal(b,1j);
        warning('tomsym:complexdata','complex data encountered. attempting to convert to tomcmplx.');
    end
    a = tomCmplx(a,spzeros(size(a)));
    y = feval('mtimes',a,b);
    return
end

% move out minus signs
if tomCmp(a,'uminus')
    y = -(-a*b);
    return
end
if tomCmp(b,'uminus')
    y = -(a*-b);
    return
end

if tomCmp(a,'setrows')
    o = operands(a);
    for i=3:2:length(o)
        o{i} = o{i}*b;
    end
    y = setrows(o{:});
elseif numel(a)==1 && numel(b)~=1
    % Scalar-matrix multiplication
    y = smtimes(a,b);
elseif numel(b)==1 && numel(a)~=1
    % Scalar-matrix multiplication
    y = smtimes(b,a);
elseif size(a,2) == size(b,1)
    % Matrix-matrix multiplication

    sz1 = size(a,1);
    sz2 = size(b,2);

    if iszero(a) || iszero(b)
        % Multiplication by zero
        y = zeros(sz1,sz2);
    elseif isnumeric(a)
        if isdiag(a) && numel(a)>1 && all(diff(getdiag(a))==0)
            y = smtimes(a(1,1),b);
        elseif size(a,1)==1 && isallone(a)
            y = sum(b,1);
        elseif nnz(a)==size(a,1) && all(sum(a,2)==1) && all(a(a~=0)==1)
            [i,j] = find(a'); %#ok
            y = submatrix(b,i,':');
        elseif size(a,1) > size(a,2) & nnz(a)==size(a,2) && all((sum(a,1)~=0)==1)
            [i,j,sc] = find(a); %#ok
            y = setrows(size(a,1),i,scalerows(sc,b));
        elseif tomCmp(b,'mtimes') && isnumeric(operand(1,b))
            % Associative multiplication of adjacent numeric matrices
            y = (a*operand(1,b))*operand(2,b);
        elseif tomCmp(b, 'vertcat')
            % Add each sub-symbol separately
            y  = spzeros(sz1,sz2);
            ix = 0;
            for i=1:nOperands(b)
                x = operand(i,b);
                y = y + a(:,ix+1:ix+size(x,1))*x;
                ix = ix+size(x,1);
            end
        elseif tomCmp(b, 'ctranspose') && tomCmpO(1,b,'horzcat')
            y  = spzeros(sz1,sz2);
            ix = 0;
            bb = operand(1,b);
            for i=1:nOperands(bb)
                x = operand(i,bb)';
                y = y + a(:,ix+1:ix+size(x,1))*x;
                ix = ix+size(x,1);
            end
        elseif isdiag(b)
            y = scalecolumns(getdiag(b),a); 
        elseif isequal(a,-1)
            y = -b;
        else
            % No simplification
            y = tomSym(mfilename,sz1,sz2,a,b);
        end
    elseif isnumeric(b)
        if isdiag(b) && numel(b)>1 && all(diff(getdiag(b))==0)
            y = smtimes(b(1,1),a);
        elseif size(b,2)==1 && isallone(b)
            y = sum(a,2);
        elseif nnz(b)==size(b,2) && all(sum(b,1)==1) && all(b(b~=0)==1)
            [i,j] = find(b); %#ok
            y = submatrix(a,':',i);
        elseif tomCmp(a,'mtimes') && isnumeric(operand(2,a))
            % Associative multiplication of adjacent numeric matrices
            y = operand(1,a)*(operand(2,a)*b);
        elseif tomCmp(a, 'horzcat')
            % Add each sub-symbol separately
            y  = spzeros(sz1,sz2);
            ix = 0;
            for i=1:nOperands(a)
                x = operand(i,a);
                y = y + x*b(ix+1:ix+size(x,2),:);
                ix = ix+size(x,2);
            end
        elseif tomCmp(a, 'ctranspose') && tomCmpO(1,a,'vertcat')
            y  = spzeros(sz1,sz2);
            ix = 0;
            aa = operand(1,a);
            for i=1:nOperands(aa)
                x = operand(i,aa)';
                y = y + x*b(ix+1:ix+size(x,2),:);
                ix = ix+size(x,2);
            end
        elseif isequal(b,-1)
            y = -a;
        elseif isdiag(a)
            y = scalerows(getdiag(a),b);
        else
            % No simplification
            y = tomSym(mfilename,sz1,sz2,a,b);
        end
    elseif isdiag(a)
        if tomCmp(a,'setdiag') && tomCmpO(1,a,'repmat') && ...
                numel(operand(1,operand(1,a)))==1
            y = smtimes(operand(1,operand(1,a)),b);
        elseif size(b,2)==1 && size(b,1)>1
            % Diagonal matrix times vector
            y = getdiag(a).*b;
        elseif tomCmp(a,'setdiag') && tomCmp(b,'setdiag')
            % Move operation inside of setdiag
            if all(size(operand(1,a)) == size(operand(1,b)))
                y = setdiag(operand(1,a).*operand(1,b));
            else
                y = setdiag(vec(operand(1,a)).*vec(operand(1,b)));
            end
        elseif isdiag(b)
            y = setdiag(getdiag(a).*getdiag(b));
        else
            y = scalerows(getdiag(a),b);
        end
    elseif isdiag(b) 
        if tomCmp(b,'setdiag') && tomCmpO(1,b,'repmat') && ...
                numel(operand(1,operand(1,b)))==1
            y = smtimes(operand(1,operand(1,b)),a);
        elseif size(a,1)==1 && size(a,2)>1
            % Vector times diagonal matrix
            y = a.*getdiag(b)';
        else
            y = scalecolumns(getdiag(b),a);            
        end
    elseif (tomCmp(a,'inv') && isequal(b,operand(1,a))) || ...
            (tomCmp(b,'inv') && isequal(a,operand(1,b)))
        % Multiplication by own inverse
        y = speye(size(a));
    elseif (tomCmp(a,'inv') && tomCmp(b,'mtimes') && isequal(operand(1,b),operand(1,a))) || ...
            (tomCmp(b,'mtimes') && tomCmpO(1,b,'inv') && isequal(operand(1,operand(1,b)),a))
        % Adjacent to own inverse
        y = operand(2,b);
    elseif (tomCmp(b,'inv') && tomCmp(a,'mtimes') && isequal(operand(2,a),operand(1,b))) || ...
            (tomCmp(a,'mtimes') && tomCmpO(2,a,'inv') && isequal(operand(1,operand(2,a)),b))
        % Adjacent to own inverse
        y = operand(1,a);
    elseif tomCmp(a,'ekron')
        args = operands(a);
        if strcmp(args{5},'I')
            args{5} = b;
        else
            args{5} = args{5}*b;
        end
        y = ekron(args{:});
    elseif tomCmp(b,'ekron')
        args = operands(b);
        if strcmp(args{4},'I')
            args{4} = a;
        else
            args{4} = a*args{4};
        end
        y = ekron(args{:});
    else
        % No simplification
        y = tomSym(mfilename,sz1,sz2,a,b);
    end
else
    error('Inner matrix dimensions must agree.');
end
