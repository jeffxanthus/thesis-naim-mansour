% function dc = sim_dc(x, Prob, varargin)
%
% TOMLAB gateway routine for the computation of the matrix of 
% constraint derivatives dc(x).
% Used when g(x) and dc(x) are computed at the same time, will call sim_gdc.m

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2003-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written May 25, 2003.   Last modified May 25, 2002.

function dc = sim_dc(x, Prob, varargin)

[g,dc] = sim_gdc(x, Prob, varargin);

% MODIFICATION LOG
%
% 030525   hkh   Written