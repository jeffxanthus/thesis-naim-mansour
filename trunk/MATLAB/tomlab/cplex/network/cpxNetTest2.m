% TOMLAB /CPLEX Network Example 2
%
% Call cplexnet with only netfile as input. 
%
% The problem is described in nexample.net

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 9.0.2$
% Written May 3, 2004.   Last modified Mar 22, 2005.

function x = cpxNetTest2

% Call with only netfile as input. Use location of cplexnetmex.mex* to
% locate the file, otherwise CPLEX will not be able to find it.

cpxnet  = which('cplexnetmex');
ps      = fileparts(cpxnet);
netfile = fullfile(ps,'nexample.net');

obj        = [];
ub         = [];
lb         = [];
tail       = [];
head       = [];
supply     = [];
callback   = [];
PriLev     = [];
BIG        = [];
cpxControl = [];
logfile    = [];
savefile   = [];
savemode   = [];

[x, slack, v, rc, f_k, Inform, Iter] = ...
    cplexnet(obj, ub, lb, tail, head, ...
    supply, callback, PriLev, BIG, ...
    cpxControl, logfile, savefile, savemode, ...
    netfile);

% MODIFICATION LOG:
%
% 040416 med  Created