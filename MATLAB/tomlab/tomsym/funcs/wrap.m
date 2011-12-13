function y = wrap(info,varargin)
% wrap - tomSym-compatible function call.
%
% y = wrap(info,...) acts as a wrapper to the function call info.fun
% making the result a tomSym symbolic object if any of the input arguments
% to FUN is. Wrap uses feval to call info.fun(...) and asserts that
% the size of y is equal to [info.sz1 info.sz2].
%
% If info.nOut is different from one, then wrap extracts the info.nOut:th
% output argument, rather than the first one.
%
% info.fun must be a text string containing the name of a function.
%
% Example: Assume that x is a tomSym object, and we want to create the
% tomSym object y = 2*f(x), where f is a scalar function which is not
% compatible with tomSym. We can then write:
%   y = 2*wrap(struct('fun','f','n',1,'sz1',1,'sz2',1),x)
%
% The function that is wrapped must return a real scalar, vector or matrix.
%
% The fields that should be present in the info structure are:
%   fun       - The name of the function to call (a character string)
%   n         - The position of the output argument that is captured
%   sz1       - The number of rows in the output
%   sz2       - The number of columns in the output
%   Jfuns     - A cell array of function names (character strings) to the
%               Jacobian matrices with respect to each input argument.
%   Jpatterns - A cell array of sparsity patterns for the Jacobians.
%   CDh       - A cell array of step lengths to use for central
%               differences.
%
% If the Jacobian matrices of the wrapped function are never requested then
% the last fields can be left out of the info structure.
%
% Jfuns can be the name of any function, or one of the following:
%   'MAD'   - Use MAD to compute derivatives
%   'FDJac' - Use one-sided finite differences with adaptive step-length
%   'CDJac' - Use central differences with fix step-length
%
% Tip: Instead of using "wrap" it is often convenient to use "feval" which
% will auto-generate the info structure for wrap if needed.
%
% See also: tomSym/feval

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-10-02 by rutquist for TOMLAB release 7.7

if ~isfield(info,'nOut')
    info.nOut = 1;
end

y = cell(1,info.nOut);
[y{:}] = feval(info.fun,varargin{:});
y = y{info.nOut};

if size(y,1)~=info.sz1 || size(y,2)~=info.sz2
    error('Wrong output size from wrapped function.');
end
