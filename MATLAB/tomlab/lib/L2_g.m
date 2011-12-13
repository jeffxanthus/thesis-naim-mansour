% L2_g.m
%
% function g=L2_g(x, Prob, varargin)
%
% L2_g computes the gradient of the rewritten L2 objective function f at x 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 13, 2002. Last modified Apr 13 2002.

function g=L2_g(x, Prob, varargin)

m = Prob.L2.m;
n = length(x)-m;
g = [zeros(n,1);x(n+1:n+m);];