function p = quickop(op, a, b)
% quickop - tomSym generic operator application for unchanged size,
% element-wise operators.
%
% c = quickop(op,a) creates an object c of the same size as a.
% c = quickop(op,a,b) checks that a and b has the same size.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009-2010 by Tomlab Optimization Inc.
% Last modified 2010-11-17 by rutquist for TOMLAB release 7.7

if(nargin==2)
    % unary operator.
    p = a;
    p.X = feval(op,a.X);

elseif(nargin==3)
    % binary operator    
    if isa(a,'tomArrayIdx')
        a = tomArray(a);
    end
    if isa(b,'tomArrayIdx')
        b = tomArray(b);
    end
    
    if ~isa(a,'tomArray')
        if numel(a)==1
            p = b;
            p.X = feval(op,a,b.X);
            return
        end
        a = tomArray(a);
        a.ni = b.ni;
    end
    if ~isa(b,'tomArray')
        if numel(b)==1
            p = a;
            p.X = feval(op,a.X,b);
            return
        end
        b = tomArray(b);
        b.ni = a.ni;
    end

    if ~isempty(a.ni) || ~isempty(b.ni)
        % Named inexes
        if isempty(a.ni) || isempty(b.ni)
            error('Either both tomArrays in a binary operation must be named, or none.');
        end

        % Expand along dimensions that are not present in both a and b
        for i=1:length(b.ni)
            if ~any(strcmp(a.ni,b.ni{i}))
                a = repmat(a,b.ni{i},b.sz(i));
            end
        end
        for i=1:length(a.ni)
            if ~any(strcmp(b.ni,a.ni{i}))
                b = repmat(b,a.ni{i},a.sz(i));
            end
        end
        
        % Make sure inexes in b are in the same order as in a.
        b = permute(b,a.ni);
        
        if ~isequal(a.sz(:),b.sz(:))
            error(['TomArray dimensions must agree. (operator ' op ')']);
        end
        
        p = a;
        p.X = feval(op,vec(a.X),vec(b.X));
        
    else
        % Unnamed indexes
        if ~isequal(a.sz(:),b.sz(:))
            error(['TomArray dimensions must agree. (operator ' op ')']);
        end
        p = a;
        p.X = feval(op,vec(a.X),vec(b.X));
    end
else
    error(['Wrong number of arguments for ' op '.']);
end
