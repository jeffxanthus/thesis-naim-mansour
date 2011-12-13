% L1_f.m
%
% function f=L1_f(x, Prob, varargin)
%
% L1_f computes the L1 objective function f in the point x

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 13, 2002. Last modified Apr 13, 2002.

function f=L1_f(x, Prob, varargin)

m = Prob.L1.m;
n = length(x)-2*m;
f = sum(x(n+1:n+2*m));