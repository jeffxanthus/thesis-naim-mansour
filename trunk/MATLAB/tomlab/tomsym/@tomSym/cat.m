function y = cat(dim,varargin)
% tomSym/cat - Overloaded function

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if numel(dim)~=1 || ~isnumeric(dim) || dim~=round(dim) || dim <= 0
    error('Dimension must be a finite integer.')
elseif dim==1
    y = vertcat(varargin{:});
elseif dim==2
    y = horzcat(varargin{:});
else
    error('High-dimensional matrices are currently not allowed with tomSym.');
end
