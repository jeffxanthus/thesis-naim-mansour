function o = vertcat(varargin)
% tomCmplx/vertcat - Concatenation of complex-valued objects

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

re = cell(size(varargin));
im = cell(size(varargin));

for i=1:length(varargin)
    re{i} = real(varargin{i});
    im{i} = imag(varargin{i});
end

o = tomCmplx(vertcat(re{:}), vertcat(im{:}));
