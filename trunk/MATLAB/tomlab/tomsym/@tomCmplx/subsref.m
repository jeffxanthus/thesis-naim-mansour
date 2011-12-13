function o = subsref(o,varargin)
% tomCmplx/subsref - Object index lookup

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-12-23 by rutquist for TOMLAB release 7.7

o.re = subsref(o.re,varargin{:});
o.im = subsref(o.im,varargin{:});
