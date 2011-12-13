% mmx_H.m
%
% function H=mmx_H(x, Prob, varargin)
%
% mmx_H computes the 2nd derivative of the minimax objective function f
% in the point x for the test problem P (Prob.P).

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 22, 1999. Last modified Apr 9, 2002.

function H=mmx_H(x, Prob, varargin)

n = length(x);
H = spalloc(n,n,0);

% MODIFICATION LOG
%
% 020409 hkh Make H sparse