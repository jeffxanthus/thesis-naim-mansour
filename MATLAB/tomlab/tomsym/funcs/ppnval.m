function y = ppnval(PPN, varargin)
% PPNVAL - N-Dimensional variation of PPVAL.
%
% Y = PPNVAL(PPN, X1, X2, ...) evaluates the piecewise polynomial given by
% PPN in the points X1, X2, ...
%
% The arguments X1, X2, ... must all be of the same size, and the returned
% Y will also have that size.
%
% If any value lies outside the range given by PPN, then polynomial
% extrapolation will be used.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-03-23 by rutquist for TOMLAB release 7.7

N = length(varargin);
if N ~= length(PPN.breaks)
    error('Dimension mismatch.');
end

b = PPN.breaks;
V = PPN.coefs;
sz = ones(1,2*N);
sz(1:length(size(V))) = size(V);
szc = sz(N+1:2*N);
sz = sz(1:N);

% Locate the interval to evaluate on
v = cell(1,N);
dx = cell(1,N);
for i=1:N
  xi = vec(full(varargin{i}));
  xi = max(xi,b{i}(1));
  xi = min(xi,b{i}(end)-eps(b{i}(end)));
  [scratch,vi] = histc(xi,b{i});
  v{i} = vi;
  dx{i} = vec(varargin{i})-vec(b{i}(vi));
end

% Evaluate polynomial
% TODO: Vectorize kk loop
V = reshape(V,prod(sz),prod(szc));
V = V(sub2ind(sz,v{:}),:);
for i=1:N
    V = reshape(V,[length(v{1}),szc(i),prod(szc(i+1:end))]);
    v1 = V(:,1,:);
    for k=2:szc(i)
        v1 = v1.*repmat(dx{i},[1,1,prod(szc(i+1:end))]) + V(:,k,:);
    end
    V = v1;
end
y = reshape(V,size(varargin{1}));
