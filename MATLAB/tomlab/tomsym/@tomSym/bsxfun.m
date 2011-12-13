function c = bsxfun(func,a,b)
% tomSym/bsxfun - Overloaded function
%
%  C = BSXFUN(FUNC,A,B) replicates the behaviour of the built-in function
%  for tomSym objects.
%
%  FUNC muts be a handle to a function that act element-wise on matrices, 
%  and which expands scalar imputs to match the size of a matrix.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(a)==1 || numel(b)==1
    c = feval(func,a,b);
else
    if (size(a,1)~=size(b,1) && size(a,1)~=1 && size(b,1)~=1) || ...
            (size(a,2)~=size(b,2) && size(a,2)~=1 && size(b,2)~=1)
        error('Matrix dimensions must agree.');
    end

    if size(a,1)==1 && size(b,1)>1
        a = repmat(a,size(b,1),1);
    elseif size(b,1)==1 && size(a,1)>1
        b = repmat(b,size(a,1),1);
    end
    if size(a,2)==1 && size(b,2)>1
        a = repmat(a,1,size(b,2));
    elseif size(b,2)==1 && size(a,2)>1
        b = repmat(b,1,size(a,2));
    end
    
    c = feval(func,a,b);
end
        
