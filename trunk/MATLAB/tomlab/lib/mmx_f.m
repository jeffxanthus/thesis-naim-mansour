% mmx_f.m
%
% function f=mmx_f(x, Prob, varargin)
%
% mmx_f computes the minimax objective function f in the point x 
% for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 22, 1999. Last modified Feb 5, 2000.

function f=mmx_f(x, Prob, varargin)
f = x(length(x));