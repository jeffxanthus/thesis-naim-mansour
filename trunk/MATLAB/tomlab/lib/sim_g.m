% function g = sim_g(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the gradient vector g(x)
% together with constraint Jacobian values dc(x), will call sim_gdc.m

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 25, 2003.   Last modified May 25, 2002.

function g = sim_g(x, Prob, varargin)

g = sim_gdc(x, Prob, varargin);

% MODIFICATION LOG
%
% 030525   hkh   Written