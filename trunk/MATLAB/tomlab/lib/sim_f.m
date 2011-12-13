% function f = sim_f(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of function values f(x)
% together with constraint values c(x), will call sim_fc.m

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 24, 2003.   Last modified May 24, 2002.

function f = sim_f(x, Prob, varargin)

f = sim_fc(x, Prob, varargin);

% MODIFICATION LOG
%
% 030524   hkh   Written