function y = mtimes(a,b)
% tomArray/mtimes - Overloaded operator
%
% Matlab's matrix multiplication operator * is, when used on tomArrays,
% only used for multiplying an array by a scalar.
%
% In order to achieve the equivalent of matrix multiplication A*B, an
% explicit sum over the element-wise product must be used, as in:
%
%   sum(A('i','j').*B('j','k'),'j')

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~((~isa(a,'tomArray') && numel(a)==1) || ...
    (~isa(b,'tomArray') && numel(b)==1))
    error(['Mtimes on tomArrays is only allowed for multiplying by scalars.' ...
        sprintf('\n') 'Use .* and sum instead.']);
end

y = quickop(mfilename,a,b);
