function [code, data, header, subfuns] = mcode(o,varargin)
% mcode - Generate m-code from a tomSym object.
%
% [code, tempD] = mcode(o) returns m-code representing the object o.
%
% Any data needed will be stored in tempD, which is referenced by the code.
%
% See also: mfile

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2011 by Tomlab Optimization Inc.
% Last modified 2011-07-14 by rutquist for TOMLAB release 7.7

if nargout <= 3
    [code, data, header] = mcode(tomSym(o),varargin{:});    
else
    [code, data, header, subfuns] = mcode(tomSym(o),varargin{:});
end
