% mmx_g.m
%
% function g=mmx_g(x, Prob, varargin)
%
% mmx_g computes the gradient of the minimax objective function f
% in the point x for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 22, 1999. Last modified Feb 5, 2000.

function g=mmx_g(x, Prob, varargin)
g = [zeros(length(x)-1,1);1];