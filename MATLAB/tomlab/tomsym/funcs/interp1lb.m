function constraints = interp1lb(x, y, xi, yi, flag)
% interp1lb - interp1 as lower bound
%
% constraints = interp1lb(x,y,xi,yi) creates a set of linear constraints
% that are equivalent to yi >= interp1(x,y,xi,'linear').
%
% constraints = interp1lb(x,y,xi,yi,'vector') is a vectorized operation
% equivalent to the loop:
% for i=1:length(xi)
%   yi(i) >= interp1(x(:,i),y(:,i),xi(i))
% end
%
% Inputs:
%
%   x   numeric vector of x-values used for interpolation
%   y   numeric vector or matrix of y-values used for interpolation
%  xi   symbol (scalar or vector) used as independent variable
%  yi   symbol (same size as xi) which will be constrained to lie on or
%       above the linear interpolant at xi.
%
% Output:
%
%  constraints - A cell array of symbolic constraints.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

if nargin<5
    flag = '';
end

if ~isnumeric(x) && isnumeric(y)
    error('x and y must be numeric');
end

if numel(y) == length(y) && numel(x) == length(x)
    mode = 1;
    
    if numel(x)~=length(x)
        error('x must be a vector');
    end
    
    if length(y) ~= length(x)
        error('x and y must have the same length');
    end
    
    x = vec(x);
    y = vec(y);
else
    if ~isempty(flag) && (flag(1)=='v' || flag(1)=='V')
       mode = 2;
       
       if length(xi)~=numel(xi)
           error('xi must be a vector.');
       end
       
       if length(x)==numel(x)
           x = repmat(x(:),1,length(xi));
       end
       
       if length(y)==numel(y)
           y = repmat(y(:),1,length(xi));
       end
       
       if size(x,2) ~= length(xi) || size(y,2) ~= length(xi)
           error('Number of columns in x and y must match the length of xi in vector mode');           
       end
       
       if size(x,1) ~= size(y,1)
           error('x and y must have the same number of rows in vector mode');
       end
       
       if ~all(size(xi)==size(yi))
           error('xi and yi must have the same size');
       end
       
    else
        mode = 3;
        error('Did you mean interp1lb(...,''vector'')?');
        % TODO
    end
end

xi = vec(xi);
yi = vec(yi);


if length(x) < 2
    error('At least two points must be specified.');
end

nx = zeros(size(x,2),1);
isconvex = true(size(x,2),1);
for k=1:size(x,2)
    xidx = find(~isnan(y(:,k)));
    xx = x(xidx,k);
    yy = y(xidx,k);
    
    if length(xx) < 2
        error('Interpolation table must contain at least two points');
    end
    
    if any(diff(xx)<=0)
        error('x must consist of strictly increasing values');
    end

    if length(xx)>2
        dy = diff(yy)./diff(xx);
        ddy = diff(dy);
        
        isused = [true; abs(ddy) >= 32*eps(dy(2:end,:)); true];
        isconvex(k) = all(ddy >= -32*eps(dy(2:end)));
    else
        isused = true(2,1);
    end
    
    nx(k) = sum(isused);
    x(1:nx(k),k) = xx(isused);
    y(1:nx(k),k) = yy(isused);

end

constraints = {};

for nn = 1:size(x,1)        
    ixc = find(isconvex & nx==nn);
    ixnc = find(~isconvex & nx==nn);
    
    if ~isempty(ixnc)
        if mode == 1
            lam = tom([],nn,length(xi));
            xs = lam'*x(1:nn);
            ys = lam'*y(1:nn);
            constraints = [constraints, {  ...
                xs == xi, ...
                ys <= yi, ...
                sos2(lam), ...
                sum(lam,1) == 1 }];
        else % mode == 2
            lam = tom([],nn,length(ixnc));
            xs = sum(x(1:nn,ixnc).*lam,1)';
            ys = sum(y(1:nn,ixnc).*lam,1)';
            constraints = [constraints, {  ...
                xs == xi(ixnc), ...
                ys <= yi(ixnc), ...
                sos2(lam), ...
                sum(lam,1) == 1 }];
        end
    end
    
    if ~isempty(ixc)
        c = cell(1,nn-1);
        if mode == 1
            for i=1:length(c)
                c{i} = ( yi >= y(i) + (xi-x(i))*(y(i+1)-y(i))/(x(i+1)-x(i)) );
            end
        else
            for i=1:length(c)
                c{i} = ( yi(ixc)' >= y(i,ixc) + (xi(ixc)'-x(i,ixc)).*(y(i+1,ixc)-y(i,ixc))./(x(i+1,ixc)-x(i,ixc)) );
            end
        end
        constraints = [constraints, c];
    end
end
