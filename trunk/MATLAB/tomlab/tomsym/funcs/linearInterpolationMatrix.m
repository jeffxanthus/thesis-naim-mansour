% linearInterpolationMatrix - Value matrix for interpPoly
% 
% M = linearInterpolationMatrix(R,X) creates a matrix M of the values of the
% the interpolating piecewise linear function defined by R, at the points X.
% Each column of M corresponds to a component of R, and each row corresponds
% to a component of X.
%
% If any value of X is very close (on the scale of machine precision) to a
% value of R, then it is modified to yield an exact match. To suppress this
% modification, use: linearInterpolationMatrix(R,X,true)

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

function M = linearInterpolationMatrix(r,x,exact)

r = full(r(:));
x = full(x(:));

% Fix numeric errors in x
if nargin<3 || ~exact
    xi = interp1(r,1:length(r),x,'nearest','extrap');
    idx = find(abs(r(xi)-x)<=4*eps*abs(x));
    x(idx) = r(xi(idx));
end

[r, order] = sort(r);
[ignore,k] = histc(x,r); %#ok
idx = (x<=r(1));
k(idx) = 1;
idx = (x>=r(end));
k(idx) = length(r)-1;
xd = (x-r(k))./(r(k+1)-r(k));
ii = [(1:length(x))'; (1:length(x))'];
jj = [k;k+1];
vv = [1-xd; xd];
idx = (vv~=0);
M = sparse(ii(idx),order(jj(idx)),vv(idx),length(x),length(r));
