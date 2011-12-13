function dpp = ppnderivative(pp,dim)
% ppderivative - The derivative of a piecewise N-D polynomial.
%
% dpp = ppnderivative(pp,dim) returns a piecewise polynomial which is the
% derivative of pp along the dimension dim.
%
% dpp = ppderivative(pp,dim, n) returns the nth derivative. Unlike the
% 1-Dimensional ppderivative, ppnderivative does not compute
% antiderivatives.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% $Id $

dpp = pp; % Starting pont.
N = length(pp.breaks);
sz = ones(1,2*N);
sz(1:length(size(pp.coefs))) = size(pp.coefs);
szc = sz(N+1:2*N);
sz = sz(1:N);

if dim<1 || dim>N
    error('Illegal dimension for derivative');
end

% Derivative of polynomial coefficients
if szc(dim)<=1
    % Derivative of piecewise constant is zero.
    dpp.coefs = zeros([sz ones(size(szc))]);
else
    % Bring the dimension to take derivative of first
    perm = [N+dim, 1:N+dim-1, N+dim+1:2*N];
    dpp.coefs = permute(dpp.coefs,perm);

    % Take derivative of polynomial
    A = spdiags((szc(dim)-1:-1:1)',0,szc(dim)-1,szc(dim));
    dpp.coefs = A*dpp.coefs(:,:);

    % Restore the shape of coefficients
    dpp.coefs = ipermute(reshape(dpp.coefs,[szc(dim)-1, sz, szc(1:dim-1) szc(dim+1:end)]),perm);
end
