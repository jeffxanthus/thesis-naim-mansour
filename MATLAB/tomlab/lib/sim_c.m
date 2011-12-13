% function c = sim_c(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the vector of constraints c(x).
% Used when f(x) and c(x) are computed at the same time, will call sim_fc.m

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2005 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 24, 2003.   Last modified May 24, 2002.

function c = sim_c(x, Prob, varargin)

[f,c] = sim_fc(x, Prob, varargin);

% MODIFICATION LOG
%
% 030524   hkh   Written