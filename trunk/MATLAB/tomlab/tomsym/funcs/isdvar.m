function y = isdvar(x)
% isdvar - Determine if a tomSym is a decision variable
%
% y = isdvar(x) returns "true" if x contains only symbols that translate
% directly to decision variables. That is: no mathematical expressions are
% allowed. (The only allowed operations are things like transpose, index
% lookup and reshape.)
%
% If x is not a tomSym, then isdvar(x) will return "false".

% Per Rutquist, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2010-2011 by Tomlab Optimization Inc.
% Last modified 2011-05-09 by rutquist for TOMLAB release 7.7

y = false;
