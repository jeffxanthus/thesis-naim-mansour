% bmi_prob:
%
% Defines Linear Semi-Definite Programming problems
% with Bilinear Matrix Inequality constraints (BMI)
%
% function [probList, Prob] = bmi_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 2004.   Last modified Jun 3, 2008.

function [probList, Prob] = bmi_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'BMI.ps Example 3',...
    'YALMIP BMI example'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_0 = []; F = []; f_opt = [];
x_L = []; x_U = []; f_Low = [];
x_min = []; x_max = [];

if P == 1
    Name='bmi.ps example 3';
    A   = [];
    b_U = [];
    b_L = [];
    c   = [0 0 1];
    % One matrix constraint, set linear part first
    SDP = [];
    % The constant matrix is stored as
    % SDP(i).Q{j} when SDP(i).Qidx(j) == 0
    SDP(1).Q{1}  = [-10 -0.5 -2 ;-0.5 4.5 0 ;-2 0 0 ];
    SDP(1).Q{2} = [ 9 0.5 0       ;  0.5  0  -3 ;  0  -3 -1 ];
    SDP(1).Q{3} = [-1.8 -0.1 -0.4 ; -0.1 1.2 -1 ; -0.4 -1 0 ];
    % Sparse is fine, too. Eventually, all the matrices are
    % converted to sparse format.
    SDP(1).Q{4} = -speye(3);
    SDP(1).Qidx = [0; 1; 2; 3];
    % Now bilinear part
    % K_12 of constraint 1 (of 1) is nonzero, so set in SDP(i).K{1}.
    SDP(1).K{1} = [0 0 2 ; 0 -5.5 3 ; 2 3 0 ];
    SDP(1).Kidx = [1 2];
    n   = length(c);
    x_L = [-5; -3; -Inf];
    x_U = [ 2;  7;  Inf];
    f_Low = [];
    % Plot ranges
    x_min=zeros(n,1);
    x_max=3*ones(n,1);
elseif P == 2
    Name='YALMIP yalmpitest.m BMI example';
    A   = [];
    b_U = [];
    b_L = [];
    c =  [0 0 0 -1];
    SDP(1).Qidx = [0; 1; 2; 3];
    SDP(1).Q{1} = [1 0; 0 1];
    SDP(1).Q{2} = [-1 0; 0 0];
    SDP(1).Q{3} = [0 -1; 0  0];
    SDP(1).Q{4} = [0  0; 0 -1];
    SDP(2).Qidx = [1; 2; 3];
    SDP(2).Q{1} = [-2 2; 0 0];
    SDP(2).Q{2} = [-6 -5; 0  4];
    SDP(2).Q{3} = [0 -3; 0 -8];
    SDP(2).Kidx = [1 4; 2 4; 3 4];
    SDP(2).K{1} = [2 0; 0 0];
    SDP(2).K{2} = [0 2; 0 0];
    SDP(2).K{3} = [0 0; 0 2];
    f_opt = -2.5;
else
    error('bmi_prob: Illegal problem number');
end

n = max( [length(c),size(A,2),length(x_L),length(x_U)] );
if isempty(x_0)
    x_0 = zeros(n,1);
end

Prob = bmiAssign(F, c, SDP, A, b_L, b_U, x_L, x_U, x_0, Name, f_Low);
Prob.f_opt = f_opt;
Prob.x_min = x_min;
Prob.x_max = x_max;
Prob.P = P;

% MODIFICATION LOG:
%
% 030120 ango Written
% 030123 ango Corrected comments
% 030124 ango Add check for x_0
% 030127 ango Change LMI(i).Q0 --> LMI(i,1).Q0
% 041117 med  xxx_prob removed and code added
% 041214 frhe Problem converted to the new LMI/BMI format
% 041229 frhe Added problem 2
% 080603 med  Switched to bmiAssign, cleaned