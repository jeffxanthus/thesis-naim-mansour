function y = numel(a,varargin) %#ok
% tomArray/numel
% Returns prod(size(a)) for tomSym a.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if nargin>2
    % TODO
    error('Subscripted numel or assignment not supported yet by tomSym.')
end

y = prod(a.sz);
