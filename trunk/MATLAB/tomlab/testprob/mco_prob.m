% mco_prob: Defines
%
% Multi-Criterium unconstrained and constrained nonlinear problems
%
% function [probList, Prob] = mco_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function [probList, Prob] = mco_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'MCO-TP1'...
    ,'MCO-TP2'...
    ,'MCO-TP3'...
    ,'MCO-TP4'...
    ,'MCO-TP5'...
    ,'MCO-TP6'...
    ,'MCO-TP7'...
    ,'MCO-TP8'...
    ,'MCO-TP9'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

JacPattern = []; t = []; weightType = []; weightY = [];
SepAlg = []; f_Low = []; A = []; b_L = []; b_U = [];
ConsPattern = []; x_min = []; x_max = []; uP = [];
c_L = []; c_U = [];

if P==1
    Name='MCO-TP1';
    % Multi Criterium optimization Test problem 1 Coello proposed by
    % Kalyanmoy Deb, P. 131-132 MC0.
    x_0=[1 0.1]';
    x_opt = [0.100000 0.600011]'; % Estimated from NPJOB
    f_opt = -27.8999999939535360; % Estimated from NPJOB
    x_L = [0.1 0.1]';
    x_U = [1.0 1.0]';
    y = zeros(2,1);
elseif P==2
    Name='MC0-TP2';
    % Multi Criterium optimization Test problem 2 Coello proposed by
    % Schaffer, P. 133 MCO
    x_0 = -4;
    x_opt = 10; % Estimated from NLPJOB
    f_opt = 0.000000000000; % Estimated from NLPJOB
    x_L = -5;
    x_U = 10;
    y = zeros(2,1);
elseif P==3
    % Multi Criterium optimization Test problem 3 Coello proposed by
    % Deb, P.134 MC0
    Name='MC0-TP3';
    x_0 = [0.1 0.1]';
    x_opt = [0.249793 0.000010]';
    f_opt = 0.750000000571891420;% Estimated from NLPJOB
    x_L = [0 -30]';
    x_U = [1 30]';
    y = zeros(2,1);
elseif P==4
    Name='MC0-TP4';
    % Multi Criterium optimization Test problem 3.1 Jonathan Writh
    % P.259 MC0
    x_0 = [-1.0 0.5]';
    x_opt = [0.981660 0.214485];
    f_opt = 1.466768566118575000;% Estimated from NLPJOB
    x_L = [-2 -2]';
    x_U = [2 2]';
    c_L = [];
    c_U = [-0.4;0.8];
    y = zeros(3,1);
elseif P==5
    Name='MC0-TP5';
    % Multi-objective Rastrigin's Problem by Zitzler & Deb, P.73 MCO
    x_0 = [0.5 0 0 0 0 0 0 0 0 0]';
    x_opt = [0 0 0 0 0 0 0 0 0 0]';
    f_opt = 0.750000135040888690;% Estimated from NLPJOB
    x_L = [0 -5 -5 -5 -5 -5 -5 -5 -5 -5]';
    x_U = [1  5  5  5  5  5  5  5  5  5]';
    y = zeros(2,1);
elseif P==6
    Name='MCO-TP6';
    % Multi-objective Griewangk's Problem presented by Kalyanmoy Deb, P.77
    x_0 = [0.34 0 0 0 0 0 0 0 0 0]';
    x_opt = [0.347810 0 0 0 0 0 0 0 0 0]';
    f_opt = 0.7500000090395866;% Estimated from NLPJOB
    x_L = [0 -512 -512 -512 -512 -512 -512 -512 -512 -512]';
    x_U = [1 511 511 511 511 511 511 511 511 511]';
    y = zeros(2,1);
elseif P==7
    Name='MCO-TP7';
    % Multi-objective Problem by Zitzler & Deb, P.78
    x_0 = [1 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5]';
    x_opt = [0 0 0 0 0 0 0 0 0 0]';
    f_opt = 1.0000000000000000;% Estimated from NLPJOB
    x_L = [0 -5 -5 -5 -5 -5 -5 -5 -5 -5]';
    x_U = [1 5 5 5 5 5 5 5 5 5]';
    y = zeros(2,1);
elseif P==8
    Name='MCO-TP8';
    % Multi-objective Problem by Srinivas and Deb, P.286 MC
    % Conbined nolinear and linear constrarined problem
    x_0 = [-5 5]';
    x_opt = [-2.5 3]';%
    f_opt = -0.25000000000000;% Estimated from NLPJOB
    x_L = [-20 -20]';
    x_U = [20 20]';
    A = [1 -3];
    b_L = [];
    b_U = -10;
    c_L = [];
    c_U = 225;
    y = zeros(2,1);
elseif P==9
    Name='MCO-TP9';
    % Multi-objective Problem Osyczka and Kundu, presented by Deb, P.288 MC
    %Conbined nolinear and linear constrarined problem
    x_0 = [3 2 2 0 1 3]';
    x_opt = [5 1 5 0 1 0]';
    f_opt = -214.0000000000;% Estimated from NLPJOB
    x_L = [0 0 1 0 1 0]';
    x_U = [10 10 5 6 5 10]';
    A = [1 1 0 0 0 0;-1 -1 0 0 0 0;1 -1 0 0 0 0;-1 3 0 0 0 0];
    b_L = [2 -6 -2 -2]';
    b_U = [];
    c_L = [-4 4]';
    c_U =[];
    y = zeros(2,1);
else
    error('mco_prob: Illegal problem number');
end

c  = 'mco_c';
dc = 'mco_dc';

if isempty(c_L) & isempty(c_U)
    c = [];
    dc = [];
end

Prob = clsAssign('mco_r','mco_J', JacPattern, x_L, x_U, Name, x_0, ...
    y, t, weightType, weightY, SepAlg, f_Low, ...
    A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
    x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.uP = uP;

if P == 7
    Prob.FUNCS.J = [];
end

% MODIFICATION LOG:
%
% 040511  med  Created
% 041117  med  xxx_prob removed and code added
% 050503  hkh  Corrected problem 2, probably still wrong, same with 7
% 060814  med  FUNCS used for callbacks instead
% 080603  med  Switched to clsAssign, cleaned