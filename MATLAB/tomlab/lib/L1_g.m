% L1_g.m
%
% function g=L1_g(x, Prob, varargin)
%
% L1_g computes the gradient of the L1 objective function f in the point x 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 13, 2002. Last modified Apr 13, 2002.

function g=L1_g(x, Prob, varargin)

m = Prob.L1.m;
n = length(x)-2*m;
g = [zeros(n,1);ones(2*m,1);];