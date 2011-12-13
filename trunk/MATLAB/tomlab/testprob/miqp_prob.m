% miqp_prob:
%
% Defines Mixed-Integer Quadratic Programming (MIQP) problems
%
% function [probList, Prob] = miqp_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2010 by Tomlab Optimization Inc. $Release: 7.5.0$
% Written June 1, 1999.   Last modified Apr 8, 2010.

function [probList, Prob] = miqp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Fletcher EQP pg 231' ...
    ,'Xpress-Optimizer Ref.Man. pg 166' ...
    ,'ibell3a' ...
    ,'itointqor' ...
    ,'ZARKO ZIVANOV' ...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_0 = []; VarWeight = []; x_opt = []; f_opt = [];
x_min = []; x_max = []; fIP = []; xIP = [];

if P == 1
    % EQP-problem. Fletcher page 231.
    x_opt=[];
    f_opt=[];
    IntVars   = logical([0 0 1]);
    Name='Fletcher EQP pg 231';
    F   = [2 0 0;0 2 0;0 0 2];
    A   = [1 2 -1;1 -1 1];
    b_L = [4 -2]';
    b_U = b_L;
    c   = zeros(3,1);
    x_0=[0 0 0]';
    x_L=[-10 -10 -10]';
    x_U=[10 10 10]';
    x_min=[0 0 -1]';
    x_max=[2 2 1]';
elseif P==2
    % The MIQP problem is from the Xpress-Optimizer Reference Manual, page 166
    %
    % min -6x(1) + 2x(1)^2 - 2x(1)x(2) + 2x(2)^2
    %
    % subject to
    %            x(1) + x(2) <= 1.9
    %            x(1),  x(2) >= 0,   x(1) integer
    c    = [-6 0]';
    Name = 'XP Ref Manual MIQP';
    F    = [4 -2;-2 4];
    A    = [1 1];
    b_L  = -Inf;
    b_U  = 1.9;
    x_L  = [0 0]';
    x_U  = [Inf Inf]';
    IntVars = [1 0];
    x_opt = [1 0.5];
    f_opt = -4.5;
elseif P==3
    Name = 'ibell3a';
    load miqp_probmat P1;
    F    = P1.F;
    c    = P1.c;
    A    = P1.A;
    b_L  = P1.b_L;
    b_U  = P1.b_U;
    x_L  = P1.x_L;
    x_U  = P1.x_U;
    IntVars = P1.IntVars;
elseif P==4
    Name = 'itointqor';
    load miqp_probmat P2;
    F    = P2.F;
    c    = P2.c;
    A    = P2.A;
    b_L  = P2.b_L;
    b_U  = P2.b_U;
    x_L  = P2.x_L;
    x_U  = P2.x_U;
    IntVars = P2.IntVars;
elseif P==5
    Name = 'ZARKO ZIVANOV';
    load miqp_probmat P3;
    F    = P3.F;
    c    = P3.c;
    A    = P3.A;
    b_L  = P3.b_L;
    b_U  = P3.b_U;
    x_L  = P3.x_L;
    x_U  = P3.x_U;
    IntVars = P3.IntVars;
else
    error('miqp_prob: Illegal problem number');
end

Prob = miqpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, ...
    IntVars, VarWeight, fIP, xIP, Name, [], [], ...
    x_min, x_max, f_opt, x_opt);

Prob.P = P;

% MODIFICATION LOG:
%
% 040102  hkh  First version, from qp_prob
% 040111  hkh  Define two simple problems from cplex and xpress distribution
% 040306  hkh  Only define variable 3 as integer in problem 1
% 041110  med  Added problem 3 and 4
% 041117  med  xxx_prob removed and code added
% 050707  med  Problem 5 added, from Zarko Zivanov
% 080603  med  Switched to conAssign, cleaned
% 100708  ango Add 5th problem to list at beginning of file
