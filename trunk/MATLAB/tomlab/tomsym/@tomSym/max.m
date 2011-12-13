function y = max(a,b,dim)
% tomSym/max - Overloaded

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-05 by rutquist for TOMLAB release 7.7

if nargin==2
    if isequal(a,b)
        y = a;
    else
        y = quickop(mfilename,a,b);
    end
else
    if nargin~=3
        if size(a,1)==1
            dim = 2;
        else
            dim = 1;
        end
    end
    if nargin>2 && ~isempty(b)
        error('Illegal syntax. Use ether max(a,b) or max(a,[],dim).');
    end
    
    if dim==1
        sz = [1, size(a,2)];
    else
        sz = [size(a,1), 1];
    end

    if size(a,dim)==1
        % Max of exactly one element
        y = a;
    elseif all(sz==[1 1]) && tomCmp(a,'max')
        y = tomSym('max',1,1,vec(operand(1,a)),[],1);
    else
        y = tomSym(mfilename,sz(1),sz(2),a,[],dim);
    end
end
