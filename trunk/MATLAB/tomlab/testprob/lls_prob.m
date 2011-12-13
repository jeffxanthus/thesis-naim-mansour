% lls_prob: Defines constrained linear least squares problems
%
% function [probList, Prob] = lls_prob(P);
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

function [probList, Prob] = lls_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'LSSOL test example'...
    ,'LSQLIN Problem Case'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

t = []; weightType = []; weightY = []; x_0 = [];
A = []; b_L = []; b_U = []; x_min = []; x_max = []; f_opt = [];
x_opt = [];

if P==1
    % Test of LSSOL. Least Squares example
    Name='LSSOL test example';
    n = 9;
    x_L = [-2 -2 -inf, -2*ones(1,6)]';
    x_U = 2*ones(9,1);
    A   = [ ones(1,8) 4; 1:4,-2,1 1 1 1; 1 -1 1 -1, ones(1,5)];
    b_L = [2    -inf -4]';
    b_U = [inf    -2 -2]';
    y = ones(10,1);
    C = [ ones(1,n); 1 2 1 1 1 1 2 0 0; 1 1 3 1 1 1 -1 -1 -3; ...
        1 1 1 4 1 1 1 1 1;1 1 1 3 1 1 1 1 1;1 1 2 1 1 0 0 0 -1; ...
        1 1 1 1 0 1 1 1 1;1 1 1 0 1 1 1 1 1;1 1 0 1 1 1 2 2 3; ...
        1 0 1 1 1 1 0 2 2];
    x_0 = 1./(1:n)';
    x_min = -ones(n,1);
    x_max =  ones(n,1);
    % x_opt estimated from LSSOL
    x_opt = [2 1.57195927 -1.44540327 -0.03700275 0.54668583 0.17512363 ...
        -1.65670447 -0.39474418  0.31002899];
    f_opt = 0.1390587318; % Estimated from LSSOL
elseif P==2
    % Matlab Newsgroup test case
    Name = 'LSQLIN Problem Case';
    load lls_probmat P1;
    y = P1.LS.y;
    C = P1.LS.C;
    x_L = P1.x_L;
    x_U = P1.x_U;
else
    error('lls_prob: Illegal problem number');
end

Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
    t, weightType, weightY, A, b_L, b_U, ...
    x_min, x_max, f_opt, x_opt);

% MODIFICATION LOG:
%
% 001106  hkh  Written
% 041117  med  xxx_prob removed and code added
% 050711  med  Added Problem 2
% 060130  med  Name added for problem 2
% 080603  med  Switched to llsAssign, cleaned