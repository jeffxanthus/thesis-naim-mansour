function o = sum(o,dim)
% tomArray/sum - Overloaded function
%
%   S = sum(X) is the sum over the first non-singelton dimension of
%   tomArray X.
%
%   S = sum(X,DIM) sums along the dimension DIM.
%
% The DIM argument can either be an integer (as in normal Matlab array
% sums) or a named dimension, for example SUM(X('i','j'),'j')

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin<2
    dim = find(o.sz>1,1,'first');
    if isempty(dim)
        o = squeeze(o);
        return
    end
end

if ischar(dim) || isa(dim,'tomArrayIdx')
    dim = find(strcmp(o.ni,char(dim)));
    if isempty(dim)
        error('The named summation index was not found in the array.');
    end
end

if o.sz(dim)==1
    % Sum over just one element.
    return
end

if length(o.sz)>1
    idxM = getIdxM(o.sz);
    n = [dim 1:dim-1 dim+1:length(o.sz)];
    idxM = reshape(permute(idxM, n),o.sz(dim),prod(o.sz)/o.sz(dim));
    o.X = sum(o.X(idxM),1);
else
    o.X = sum(o.X(:));
end

o.sz(dim) = [];
if ~isempty(o.ni)
    o.ni = {o.ni{1:dim-1}, o.ni{dim+1:end}};
end

    
