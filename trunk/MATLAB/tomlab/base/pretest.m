% pretest.m
%
% function [A, b, conType, x_L, x_U] = pretest(P)
%
% pretest.m is a test procedure for preSolve.m .
%
% INPUT PARAMETERS
% P         Problem number, default P=1. Problems 1-9 defined
%
% OUTPUT PARAMETERS
%   A       Linear constraint matrix
%   b_L     Lower bounds for linear constraints
%   b_U     Upper bounds for linear constraints
%   x_L     Lower bounds for x
%   x_U     Upper bounds for x

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written Jan 25, 1998.  Last modified Jan 25, 2004.

function [A, b_L, b_U, x_L, x_U] = pretest(P)

if nargin < 1
   P = 1;
end

x_L = [];
x_U = [];
b_L = [];
b_U = [];
A   = [];

if P == 1
   a1 = [ 1  1  1  1 ];
   a2 = [ 0  1  1  1 ];
   a3 = [ 0  0  1  1 ];
   a4 = [ 0  0  0  1 ];
   A  = [a1;a2;a3;a4];
   b_L  = [1 2 3 4];
   b_U  = [1 2 3 4];
   %conType = [11 11 11 11];
   %x_L = [-100 -100 -100 -100];
   %x_U = [ 100  100  100  100];
elseif P == 2
   a1 = [ 0  0  1e-2];
   a2 = [ 1  2    3  ];
   A  = [a1;a2];
   b_L = [0.5 4];
   b_U = [0.5 4];
   %conType = [11 11];
   x_L = -[100 100 100];
   x_U =  [100 100 100];
elseif P == 3
   a1 = [ 1  1 ];
   a2 = [ 1 -1 ];
   a3 = [-1  1 ];
   a4 = [-1 -1 ];
   A  = [a1;a2;a3;a4];
   b_U = [3 3 3 3];
   %conType = [13 13 13 13];
   x_L = [-1 -inf ];
   x_U = [ 500  inf ];   
elseif P == 4
   a1 = [ 0  1  2  3 ];
   a2 = [ 2  0  0  3 ];
   a3 = [ 1  2  3  4 ];
   a4 = [ 0  0  2  0 ];
   A  = [a1;a2;a3;a4];
   b_L  = [ 1  2  3  4 ];
   b_U  = [ 1  2  3  4 ];
   %conType = [11 11 11 11];
   x_L = [-100 -100 -100 -100 ];
   x_U = [ 100  100  100  100 ];
A=[A;A];
b_L=[b_L b_L];
b_U=[b_U b_U];
%conType=[conType conType];   

elseif P == 5
   A = [  1   1   1   1   1   1   1	
	      15   4   2   4   2   1   3
  	       3   5   8   2   6   1   0
	       2   4   1   2   2   1   0
	       2   3   0   0   1   0   0
	      70  75  80  75  80  97   0
	       2   6   8  12   5   1  97
	       2   6   8  12   5   1  97
	       1   0   0   0   0   0   0
	       0   1   0   0   0   0   0
	       0   0   1   0   0   0   0
	       0   0   0   1   0   0   0
          0   0   0   0   1   0   0 ];
   b_L     = [15 -inf -inf -inf -inf  1105  206 -inf -inf -inf -inf -inf -inf ];
   b_U     = [15  48   66   34   30    inf  inf  256   2    25  4     6    15 ];   
  %conType = [11  13   13   13   13    12   12    13   13   13  13   13    13 ];
   x_L = zeros(1,7);
   x_U = 1e9.*ones(1,7);
elseif P == 6
   A = hilb(3);
   b_L = [1 2 3];
   b_U = [1 2 3];
  %conType = [11 11 11];
   x_L = [-1000 -1000 -1000];
   x_U = [ 1000  1000  1000];   
elseif P == 7 % mbk test 2 from cls-problems
   A = [5.0    1.0
        1.25   1.0
        0.5    1.0
        0.25   1.0];
   b_L = [0.000000   4.000000   4.000000   3.750000];
   %b_U = inf* ones(4,1);
   conType = [12 12 12 12];
   x_L = [ 1.0 -1E2 ]; 
   x_U = [ 1.5  1E2 ];
elseif P == 8
   A = [1  2  3  4
        5  6  7  8
        0 0 0 1 ];
     b_L = [2 3 4];
     b_U = [2 13 2332];
     
   x_L = []; 
   x_U = [];
elseif P == 9
   A = eye(4);
     b_L = [1 2 3 4];
     b_U = [5 6 7 inf];
     
   x_L = []; 
   x_U = [];
else
   fprintf('\n Problem %d not defined',P);
   return;
end

Prob.A = A;
Prob.b_L = b_L;
Prob.b_U = b_U;
Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.MIP = [];
Prob.PriLevOpt = 1;

Prob = preSolve(Prob);

A   = Prob.A;
b_L = Prob.b_L;
b_U = Prob.b_U;
x_L = Prob.x_L;
x_U = Prob.x_U;

% MODIFICATION LOG:
%
% 990622  hkh  Major parts written
% 040125  hkh  Revisions for v4.2. Define fields MIP and PriLevOpt
