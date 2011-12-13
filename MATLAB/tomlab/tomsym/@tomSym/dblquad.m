function y = dblquad(fun,a1,b1,a2,b2,x1,x2)
% tomSym/dblquad - Numeric double quadrature
%
% y = dblquad(fun,a1,b1,a2,b2,x1,x2) - Computes the double integral
% over the rectangle a1 < x1 < b1, a2 < x2 <b2 of fun, with 
% respect to x1 and x2. fun should be a tomSym expressoin
%
% y = dblquad(fun,a1,b1,a2,b2) - where fun is a function handle, calls that
% function with symbolic inputs to obtain a tomSym function. (This syntax
% is compatible with the Matlab function that is overloaded.)
%
% If fun contains more than one symbol, then the integral is not computed
% until all the remaining symbols have been substituted by constants.

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2009 by Tomlab Optimization Inc.
% Last modified 2009-08-21 by rutquist for TOMLAB release 7.7

if ~isa(fun,'tomSym') && nargin<=5
    x1 = tom([],1,1);
    x2 = tom([],1,1);
    fun = feval(fun,x1,x2);
end

y = quad(quad(fun,a1,b1,x1),a2,b2,x2);

