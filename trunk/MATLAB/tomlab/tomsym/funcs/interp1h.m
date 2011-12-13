function yi = interp1h(x,y,xi)
% tomSym/interp1h - Piecewise constant (hold) interpolation
%
% yi = interp1h(x,y,xi) interpolates the piecewise constant function given
% by (x,y) at the points xi. This is similar to interp1(...,'nearest'),
% except the boundaries are given by x, so for example, y(1) is returned
% for any xi in the intervall from x(1) to x(2).

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~any(size(x)==1)
    error('X must be a vector.');
end

if size(x,2)>1
    x = x';
end

if size(y,1)==1 && size(x,1)>1 && size(y,2)==size(x,1)
    y = y';
end

if size(x,1)~=size(y,1)
    error('Input arguments x and y must have the same length.');
end

[x, order] = sort(full(x(:)));
[ignore,k] = histc(full(xi(:)),x); %#ok
k(xi<=x(1)) = 1;
k(xi>=x(end)) = length(x);

yi = y(order(k),:);
if size(y,2)==1
    yi = reshape(yi,size(xi));
end
