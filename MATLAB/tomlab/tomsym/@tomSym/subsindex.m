function [varargout] = subsindex(varargin) %#ok
% tomSym/subsindex - Not possible

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-26 by rutquist for TOMLAB release 7.7

error(['Cannot use a tomSym object directly as an index (due to Matlab syntax rules). ' 10 ...
    '- Use lookup or submatrix instead. (See "help lookup")']);
