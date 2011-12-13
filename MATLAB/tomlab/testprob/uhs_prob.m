% Defines unconstrained optimization problems (with simple bounds)
%
% uhs_prob: Defines unconstrained optimization problems (with simple bounds)
%
% function [probList, Prob] = uhs_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2011 by Tomlab Optimization Inc. $Release: 7.8.0$
% Written June 1, 1999.   Last modified July 25, 2011.

function [probList, Prob] = uhs_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'HS 1'...
    ,'HS 2'...
    ,'HS 3'...
    ,'HS 4'...
    ,'HS 5'...
    ,'HS 6'...
    ,'HS 7'...
    ,'HS 8'...
    ,'HS 9'...
    ,'HS 10'...
    ,'HS 11'...
    ,'HS 12'...
    ,'HS 13'...
    ,'HS 14'...
    ,'HS 15'...
    ,'HS 16'...
    ,'HS 17'...
    ,'HS 18'...
    ,'HS 19'...
    ,'HS 20'...
    ,'HS 21'...
    ,'HS 22'...
    ,'HS 23'...
    ,'HS 24'...
    ,'HS 25'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

HessPattern = []; pSepFunc = []; f_Low = [];

if P == 1
    Name='HS 1';
    x_0 = [-2,1]';
    x_L = [-inf;-1.5];
    x_U = [];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = [-2 -2];
    x_max = [ 2  2];
elseif P == 2
    Name='HS 2';
    x_0 = [-2,1]';
    x_L = [-inf;1.5];
    x_U = [];
    a = sqrt(598/1200);
    b = 400*a^3;
    x_opt = [2*a*cos(1/3*acos(1/b));1.5];
    f_opt = 0.0504261879;
    f_Low = 0 ;
    x_min = [-2 -2];
    x_max = [ 2  2];
elseif P == 3
    Name='HS 3';
    x_0 = [10,1]';
    x_L = [-inf;0];
    x_U = [];
    x_opt = [0;0];
    f_opt = 0;
    x_min = [-20 -20];
    x_max = [ 20  20];
elseif P == 4
    Name='HS 4';
    x_0 = [1.125,0.125]';
    x_L = [1;0];
    x_U = [];
    x_opt = [1;0];
    f_opt = 8/3;
    x_min = [1 0];
    x_max = [2 1];
elseif P == 5
    Name='HS 5';
    x_0 = [0,0]';
    x_L = [-1.5;-3];
    x_U = [4;3];
    x_opt = [-pi/3+0.5;-pi/3-0.5];
    f_opt = -sqrt(3)/2-pi/3;
    x_min = x_L;
    x_max = x_U;
elseif P == 6
    Name='HS 38';
    x_0 = [-3,-1,-3,-1]';
    x_L = [-10,-10,-10,-10]';
    x_U = [10,10,10,10]';
    x_opt = [1,1,1,1]';
    f_opt = 0;
    x_min = x_L;
    x_max = x_U;
elseif P == 7
    Name='HS 45';
    x_0 = [2,2,2,2,2]';
    x_L = [0,0,0,0,0]';
    x_U = (1:5)';
    x_opt = (1:5)';
    f_opt = 1;
    x_min = x_L;
    x_max = x_U;
elseif P == 8
    Name='HS 110';
    x_0 = [9,9,9,9,9,9,9,9,9,9]';
    x_L = [2.001,2.001,2.001,2.001,2.001 ...
        2.001,2.001,2.001,2.001,2.001]';
    x_U = [9.999,9.999,9.999,9.999,9.999 ...
        9.999,9.999,9.999,9.999,9.999]';
    x_opt = [9.35025655,9.35025655,9.35025655,9.35025655,9.35025655,...
        9.35025655,9.35025655,9.35025655,9.35025655,9.35025655]';
    f_opt = -45.77846971;
    x_min = x_L;
    x_max = x_U;
elseif P == 9
    Name='HS 206';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 10
    Name='HS 207';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 11
    Name='HS 208';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 12
    Name='HS 209';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 13
    Name='HS 210';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 14
    Name='HS 211';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 15
    Name='HS 212';
    x_0 = [2,0]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [0 0];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 16
    Name='HS 213';
    x_0 = [3,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 17
    Name='HS 214';
    x_0 = [-1.2,1]';
    x_L = [-inf;-inf;];
    x_U = [inf;inf;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 18
    Name='HS 229';
    x_0 = [-1.2,1]';
    x_L = [-2;-2;];
    x_U = [2;2;];
    x_opt = [1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 19
    Name='HS 247';
    x_0 = [-1,1,1]';
    x_L = [0.1;-inf;-2.5];
    x_U = [inf;inf;7.5];
    x_opt = [1 0 0];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 20
    Name='HS 256';
    x_0 = [3,-1,0,1]';
    x_L = [-inf;-inf;-inf;-inf];
    x_U = [inf;inf;inf;inf];
    x_opt = [0 0 0 0];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 21
    Name='HS 257';
    x_0 = [-3,-1,-3,-1]';
    x_L = [0;-inf;0;-inf];
    x_U = [inf;inf;inf;inf];
    x_opt = [1 1 1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 22
    Name='HS 258';
    x_0 = [-3,-1,-3,-1]';
    x_L = [-inf;-inf;-inf;-inf];
    x_U = [inf;inf;inf;inf];
    x_opt = [1 1 1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 23
    Name='HS 259';
    x_0 = [0,0,0,0]';
    x_L = [-inf;-inf;-inf;-inf];
    x_U = [inf;inf;inf;inf];
    x_opt = [1 1 1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 24
    Name='HS 260';
    x_0 = [-3,-1,-3,-1]';
    x_L = [-inf;-inf;-inf;-inf];
    x_U = [inf;inf;inf;inf];
    x_opt = [1 1 1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
elseif P == 25
    Name='HS 273';
    x_0 = [0,0,0,0,0,0]';
    x_L = [-inf;-inf;-inf;-inf;-inf;-inf];
    x_U = [inf;inf;inf;inf;inf;inf];
    x_opt = [1 1 1 1 1 1];
    f_opt = 0;
    f_Low = 0 ;
    x_min = x_L;
    x_max = x_U;
else
    error('uc_prob: Illegal problem number')
end

% Define the Prob
Prob = conAssign('uhs_f','uhs_g','uhs_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, [], [], [], [], [],...
    [], [], [], [], x_min, x_max, f_opt, x_opt);
Prob.P  = P;

% MODIFICATION LOG:
%
% 990601  mbk  First version.
% 030211  hkh  Adapted for v4.0
% 040430  ango Fixed several inf:inf that should be inf;inf
% 041117  med  xxx_prob removed and code added
% 080603  med  Switched to conAssign, cleaned
% 110725  hkh  Wrong dimension on x_opt for problem 8
