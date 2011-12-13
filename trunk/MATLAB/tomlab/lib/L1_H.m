% L1_H.m
%
% function H=L1_H(x, Prob, varargin)
%
% L1_H computes the 2nd derivative of the L1 objective function f 
% in the point x 

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Apr 13, 2002. Last modified Apr 13, 2002.

function H=L1_H(x, Prob, varargin)

n = length(x);
H = spalloc(n,n,0);

% MODIFICATION LOG
%
% 020413 hkh Written