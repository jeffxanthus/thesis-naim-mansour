% function g = oqnlp_g(x, Prob, varargin)
%
% TOMLAB interface routine to avoid numerical gradient differentiation of
% integer variables in MINLP problems

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Nov 13, 2004.   Last modified Nov 13, 2004.

function g = oqnlp_g(x, Prob, varargin)

n          = Prob.N;
IntVars    = Prob.MIP.IntVars;
g          = NaN*ones(n,1);
g(IntVars) = 0;

% MODIFICATION LOG
%
% 041113  med  Written