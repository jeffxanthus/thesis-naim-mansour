% L2_f.m
%
% function f=L2_f(x, Prob, varargin)
%
% L2_f computes the rewritten L2 objective function f in the point x 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 13, 2002. Last modified Apr 13, 2002.

function f=L2_f(x, Prob, varargin)

m = Prob.L2.m;
n = length(x)-m;
r = x(n+1:n+m);
f = 0.5*(r'*r);