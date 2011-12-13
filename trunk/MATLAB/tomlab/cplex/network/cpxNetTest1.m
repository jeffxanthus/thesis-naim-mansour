% TOMLAB /CPLEX Network Example 1
%
% Example: Network Optimizer in the Interactive Optimizer:
% ILOG CPLEX Manual.
%
% See the TOMLAB /CPLEX Manual for a visual problem description.
% This example is based on a network where the aim is to minimize
% cost and where the flow through the network has both cost and capacity.
% The nodes (8) and arcs (14) are listed below. The lower/upper bound and
% cost for each arc is listed below. N1 and N5 are sources, N4 and N8 are
% sinks.
%
% Nodes:
% N1: Supply + 20, 1 arcs.
% N2: Transport node, 5 arcs.
% N3: Transport node, 3 arcs.
% N4: Demand -15, 5 arcs.
% N5: Supply +5, 4 arcs.
% N6: Transport node, 6 arcs.
% N7: Transport node, 2 arcs.
% N8: Demand -10, 2 arcs.
%
% Arcs (lower bound, upper bound, cost, connection):
% A1:  18,   24,  $3, N1->N2
% A2:  0,    25,  $3, N2->N3
% A3:  12,   12,  $4, N3->N4
% A4:  0,    10,  $3, N4->N7
% A5:  0,    9,   $5, N7->N6
% A6:  -inf, inf, $6, N6->N8
% A7:  10,   20,  $7, N5->N8
% A8:  0,    10,  $4, N5->N2
% A9:  0,    5,   $2, N3->N2
% A10: 0,    15,  $6, N4->N5
% A11: 0,    10,  $5, N4->N6
% A12: 0,    11,  $4, N6->N4
% A13: 0,    6,   $3, N6->N5
% A14: 0,    inf, $6, N2->N6

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2004-2005 by Tomlab Optimization Inc., $Release: 9.0.2$
% Written Apr 17, 2004.   Last modified Mar 22, 2005.

function x = cpxNetTest1

%%%%%%%%%%%%% ARC INFORMATION %%%%%%%%%%%%%%%%%%%

% The cost vector for the arcs.
obj =  [3  3  4  3  5 6    7  4  2 6  5  4  3 6]';

% The lower bounds for the arcs.
lb  =  [18 0  12 0  0 -inf 10 0  0 0  0  0  0 0]';

% The upper bounds for the arcs.
ub  =  [24 25 12 10 9 inf  20 10 5 15 10 11 6 inf]';

% Start connection or tails.
tail = [1  2  3  4  7 6    5  5  3 4  4  6  6 2]';

% End connections or heads.
head = [2  3  4  7  6 8    8  2  2 5  6  4  5 6]';

%%%%%%%%%%%%% NODE INFORMATION %%%%%%%%%%%%%%%%%%

supply = [20 0 0 -15 5 0 0 -10]';

%%%%%%%%%%%%% SOLVER CALL %%%%%%%%%%%%%%%%%%%%%%%

cpxControl = [];

[x] = cplexnet(obj, ub, lb, tail, head, supply, [], [], [], cpxControl);