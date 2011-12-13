% splineInterpolationMatrix - Value matrix for interpPoly
% 
% M = splineInterpolationMatrix(x,xi) creates a matrix M of the values of
% the the interpolating polynomials defined by the interpolation points x, 
% at the evaluation points xi.
%
% This means that yi=M*y gives yi so that (xi,yi) lie on the  the points
% that interpolate (x,y).
%
% M = splineInterpolationMatrix(x,xi,N) gives a similar matrix, but which
% works with the N:th derivative of the spline curve. (Negative N gives the
% antiderivative (integral) of the curve.)
%
% Each column of M corresponds to a specific interpolation point, 
% and each row corresponds to an evaluation point.
%
% If any value of xi is very close (on the scale of machine precision) to a
% value of x, then it is modified to yield an exact match. To suppress this
% modification, use: proptInterpolationMatrix(x,xi,true)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

function M = splineInterpolationMatrix(x,xi,n,exact)

if nargin<3
    n = 0;
end

x = full(x(:));
if isnumeric(xi)
    xi = full(xi(:));
    % Fix numeric errors in xi
    if nargin<4 || ~exact
        i = interp1(x,1:length(x),xi,'nearest','extrap');
        idx = find(abs(x(i)-xi)<=4*eps*abs(xi));
        xi(idx) = x(i(idx));
    end
else
    xi = vec(xi);
end

Mc = cell(1,length(x));

% Compute the matrix
for i = 1:length(x)
    y = zeros(size(x));
    y(i) = 1;
    pp = spline(x, y);
    pp = ppderivative(pp,n);
    Mc{i} = ppval(pp,xi);
end

M = horzcat(Mc{:});

if isnumeric(M) && nnz(M)<numel(M)
    M = sparse(M);
end
