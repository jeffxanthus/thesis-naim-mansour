% glcIP_prob.m
%
% Defines constrained global optimization problems with Integer Variables.
%
% function [probList, Prob] = glcIP_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob     Problem Structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2008-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written Sept 25, 2008.   Last modified Sept 29, 2008.

function [probList, Prob] = glcIP_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
     'Floudas-Pardalos 12.2 TP 1'...
    ,'Floudas-Pardalos 12.2 TP 2'...
    ,'Floudas-Pardalos 12.2 TP 3'...
    ,'Floudas-Pardalos 12.2 TP 4'...
    ,'Floudas-Pardalos 12.2 TP 5'...
    ,'Floudas-Pardalos 12.2 TP 6'...
    ,'Floudas-Pardalos 3.3 TP 2 IP'...
    ,'Floudas-Pardalos 3.4 TP 3'...
    ,'Floudas-Pardalos 3.5 TP 4 IP'...
    ,'Floudas-Pardalos 4.10 TP 9 IP'...
    ,'Kocis & Grossmann 1989'...
    ,'Kocis & Grossmann 1998'...
    ,'Lee & Grossmann Disjunctive B&B'...
    ,'Kesavan et al. 2004 D'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES
% NOT USED:
%   ,'Floudas-Pardalos 12.2 TP 4 GAMS version'...
%   ,'Floudas-Pardalos 12.2 TP 4 GAMS Int 1st '...
%   ,'Floudas-Pardalos 12.2 TP 1 IntV last'...


if isempty(P)
    return;
end

f_Low = [];
%x_0  = []; % Currently not used
user = [];

if P == 1   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 1';
    % Kocis and Grossmann (1988)
    x_L = [ 0 0 0  0   0  ]';
    x_U = [ 1 1 1 1.6 1.5 ]'; % Bnd x_1 from 3rd ineq., Bnd x_2 from 4th ineq
    A = [1 0 0 1 0; 0 1 0 0 1.333; -1 -1 1 0 0];
    b_L=[-inf -inf -inf]';
    b_U=[1.6 3 0];
    c_L = [1.25  3 ]';
    c_U = [1.25  3 ]';
    x_opt = [0 1 1 1.12 1.31];
    f_opt = 7.6672;
    IntVars = 1:3;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 2   %Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 2';
    % Floudas (1995) example 6.6.5
    x_L = [ 0 0.2 -2.22554 ]';
    x_U = [ 1 1   -1      ]';
    A = [ 1.1 0 1 ; -1.2 1 0 ];
    b_L=[-inf -inf]';
    b_U=[ -1   0.2];
    c_L = -inf;
    c_U = 0;
    % This is the optimal solution given, does not seem OK
    % x_opt = [1 0.94149 -2.1 ];
    % f_opt = 1.07456710;
    % Solution obtained with glcCluster/multiMin/snopt
    % To get more digits in f, must have analytic derivatives
    x_opt   = [1 0.941937 -2.1 ];
    f_opt   = 1.076543;
    IntVars = 1;
    x_min   = x_L;
    x_max   = x_U;
    x_0     = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 3   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 3';
    % Yuan et al. (1988)
    x_L = zeros(7,1);
    x_U = [1.2 1.8 2.5 1 1 1 1]'; % Added bnds x(1:3) from lin.eq 2-4
    A = [1 1 1 1 1 1 0; eye(3,3),eye(3),zeros(3,1);1 0 0 0 0 0 1];
    b_L=-inf*ones(5,1);
    b_U=[  5 1.2 1.8 2.5 1.2]';
    c_L=-inf*ones(4,1);
    c_U = [  5.5 1.64 4.25 4.64]';
    x_opt = [0.2 0.8 1.908 1 1 0 1];
    f_opt = 4.5796;
    IntVars = 4:7;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 4   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 4';
    % Berman and Ashrafi (1993)
    x_L = zeros(11,1);
    x_U = ones(11,1);
    A = [zeros(4,3),[-1 -1 -1, zeros(1,5); 0 0 0 -1 -1 -1 0 0; ...
        zeros(1,6),-1,-1; 3 1 2 3 2 1 3 2]];
    b_L = -inf*ones(4,1);
    b_U = [  -1 -1 -1 10]';
    c_L = ones(3,1);
    c_U = ones(3,1);
    x_opt = [0.97, 0.9925, 0.98, 0 1 1 1 0 1 1 0];
    f_opt = -0.9434705;
    IntVars = 4:11;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 5   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 5';
    % Porn et al. (1987)
    x_L = ones(2,1);
    x_U = 5*ones(2,1);
    A   = [-1 -2 ; -3 1;4 -3];
    b_L = -inf*ones(3,1);
    b_U = [  -5 1 11]'; % Wrong sign in book for 1st constraint
    % Bounds as given originally - gives different solution [1 1];
    c_L = -inf;
    c_U =   -24; % Wrong sign in book
    x_opt = [3 1];
    f_opt = 31;
    IntVars = 1:2;
    x_min = x_L;
    x_max = x_U;
    x_0 = x_L;
elseif P == 6   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 12.2 TP 6';
    % Porn et al. (1987)
    x_L = ones(2,1);
    x_U = [10 6]';
    A   = [1 -1 ; 3 2];
    b_L = -inf*ones(2,1);
    b_U = [3 24]';
    c_L = -inf;
    c_U =   39;
    x_opt = [4 1];
    f_opt = -17;
    IntVars = logical([0 1]); % Use logical format
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 7   % Problem from old glc_prob
    % Added integer constraints to FP 3.3 TP 2 and created an IP.
    Name = 'Floudas-Pardalos 3.3 TP 2 IP';
    x_L = [78 33 27 27 27]';
    x_U = [102 45 45 45 45]';
    b_L = []; b_U = []; A = [];
    c_L = [-85.334407;9.48751;10.699039];
    c_U = [6.665593;29.48751;15.699039];
    % x_opt = [78 33 29.9953 45 36.7758];
    %IntVars = [1 2];    % 1:5 ???
%     IntVars = 1:5;
%     x_opt = [81 33 30 45 36];
%     f_opt = -30512.4499954;
    IntVars = 3:5;
    x_opt = [80.862369 33.024691 30 45 36];
    f_opt = -30521.723332775728;
%     IntVars = 4:5;
%     x_opt = [78 33 30.021353 42 38];
%     f_opt = -30577.349859767826;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 8   % Problem from old glc_prob
    Name = 'Floudas-Pardalos 3.4 TP 3';
    tol = 0; % Avoid changing bounds
    x_L   = [0,  0,    1,  0-tol,  1,       0     ]';
    x_U   = [6,  6,  5+tol, 6,     5+tol,   10+tol]';
    A = [ 1 -3 0 0 0 0
        -1  1 0 0 0 0
        1  1 0 0 0 0];
    b_L = [-inf -inf  2 ]';
    b_U = [  2    2   6 ]';
    c_L = [4 4]';
    c_U = [];
    x_opt = [5 1 5 0 5 10];
    IntVars = 3:4;
    f_opt = -310;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 9   % Problem from old glc_prob
    % Added integer constraint to FP 3.5 TP 4 and created an IP.
    Name = 'Floudas-Pardalos 3.5 TP 4 IP';
    x_L = [0 1 0]'; % Lower bound on x(2) added by nhq to create IP-problem
    x_U = [2 2 3]'; % Upper bounds on x(2) added by hkh, from lin.eq.2
    A = [1 1 1
        0 3 1];
    b_L = -inf*ones(2,1);
    b_U = [4 6]';
    bvtmp = [3 1 2]';
    rtmp  = [1.5 -0.5 -5]';
    c_L = 0.25*bvtmp'*bvtmp-rtmp'*rtmp;
    c_U = [];
    IntVars = 2;
    x_opt = [2 1 0.5];
    f_opt = -3.5;
%     IntVars = 2:3;
%     x_opt = [1.634005 1 1];
%     f_opt = -3.268009054626948;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 10  % Problem from old glc_prob
    % Added integer constraint to FP 4.10 TP 9 and created an IP.
    Name = 'Floudas-Pardalos 4.10 TP 9 IP';
    x_L = [0 0]';
    x_U = [3 4]';
    A = []; b_L=[]; b_U=[];
    c_L = [-2 -36]';
    c_U = [];
    %IntVars = 1, 2 or [1 2] ???
    % 1 or [1,2] -> x_opt = [2,2], f_opt = -4
    % 2          -> x_opt = [2.3661 3], f_opt = -5.3661
%     IntVars = 1; % Intvars = 1:2;
%     x_opt = [2 2];
%     f_opt = -4;
    IntVars = 2;
    x_opt = [2.366094260698742 3];
    f_opt = -5.366094260698741;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
elseif P == 11  % Problem from minlp_prob
    Name = 'Kocis & Grossmann 1989';
    % Variables [v,x,c,p,y1,y2];
    IntVars = [3 4];
    x_L = [  0  0  0 0 ]';
    x_U = [ 10 20  1 1 ]';
    A   = [0 0 1 1 ];
    b_L = 1;
    b_U = 1;
    c_L = -inf*ones(4,1);
    user.M = 100; % use big-M method
    x_0 = [0 0 0 1];
    % This solution cannot be optimal, nonlinear constraint 2 is infeasible
    % Computed solution, with fixed integer
    x_opt = [3.514237,13.427995,1,0];
    f_opt = 99.239635053654084;
    c_U   = [0;0;1;0];
    x_min = [ 0  0 0 0 ];
    x_max = [10 20 1 1 ];
elseif P == 12  % Problem from minlp_prob
    Name = 'Kocis & Grossmann 1998';
    IntVars = [3 4 5];
    BIGBND = 1E4;
    % May get division by zero if not x_L(2)==small but >0
    x_L = [ 0      1/BIGBND   0 0 0 ]';
    x_U = [ BIGBND  BIGBND    1 1 1 ]';
    A = [1 0     1  0 0 ; ...
        0 1.333 0  1 0 ; ...
        0 0    -1 -1 1 ];
    b_L = -inf*ones(3,1);
    b_U = [1.6 ; 3 ; 0];
    c_L = [1.25;3];
    c_U = c_L;
    % Setting x_0 = ones(5,1); gives a non-feasible point with respect to
    % linear constraints. Adjusting this point gives (0.6,1,1,1,1)
    % However from this initial point MINLP will converge to a
    % nonglobal solution
    %
    % If disabling the feasibility test, then MINLP converges
    x_0 = ones(5,1);
    x_opt = [1.12,1.31,0,1,1];
    f_opt = 7.6672;
    x_min = [-1 -1 0 0 0];
    x_max = [ 1  1 1 1 1];
elseif P == 13  % Problem from minlp_prob
    Name = 'Lee & Grossmann Disjunctive B&B';
    % Variables [x(1),x(2),y1,y2,y3];
    IntVars = [3 4 5];
    x_L = [  0  0  0 0 0 ]';
    x_U = [  8  8  1 1 1 ]';
    A   = [0 0 1 1 1 ];
    b_L = 1;
    b_U = 1;
    c_L = -inf*ones(3,1);
    c_U = [0;0;0];
    user.M = 30; % use big-M method
    x_0 = [0 0 0 0 1];
    x_opt = [3.293,1.707, 0, 1, 0];
    f_opt = 1.172;
    x_min = x_L;
    x_max = x_U;
elseif P == 14  % Problem from minlp_prob
    Name = 'Kesavan et al. 2004 D';
    IntVars = (1:3)';
    x_L = [0 0 0  1 1]';
    x_U = [1 1 1 10 6]';
    A   = [0 0 0 1 -1;...
        0 0 0 3  2;...
        1 2 4 0 -1];
    b_L = [-inf; -inf; 0];
    b_U = [3; 24; 0];
    c_L = -inf*ones(1,1);
    c_U = 39;
    x_0 = [0 0 0 1 1];
    x_opt = [1 0 0 4 1];
    f_opt = -17;
    x_min = x_L;
    x_max = x_U;
elseif P == 416 % Problem from old glc_prob - NOT USED
    Name = 'Floudas-Pardalos 12.2 TP 4 GAMS version';
    % Berman and Ashrafi (1993)
    % GAMS Version has changed upper bounds on x(1:3) compared to P=22
    % The constraints are log transformed
    x_L = zeros(11,1);
    x_U = ones(11,1);
    x_U(1:3) = [0.997;0.9985;0.9988];
    A = [zeros(4,3),[-1 -1 -1, zeros(1,5); 0 0 0 -1 -1 -1 0 0; ...
        zeros(1,6),-1,-1; 3 1 2 3 2 1 3 2]];
    b_L = -inf*ones(4,1);
    b_U = [  -1 -1 -1 10]';
    c_L = zeros(3,1);
    c_U = zeros(3,1);
    x_opt = [0.97, 0.9925, 0.98, 0 1 1 1 0 1 1 0];
    f_opt = -0.9434705;
    IntVars = 4:11;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 417 % Problem from old glc_prob - NOT USED
    Name = 'Floudas-Pardalos 12.2 TP 4 GAMS Int 1st ';
    % Berman and Ashrafi (1993)
    % The constraints are log transformed
    x_L = zeros(11,1);
    x_U = ones(11,1);
    x_U(9:11) = [0.997;0.9985;0.9988];
    A = [[-1 -1 -1, zeros(1,5); 0 0 0 -1 -1 -1 0 0; ...
        zeros(1,6),-1,-1; 3 1 2 3 2 1 3 2],zeros(4,3)];
    b_L = -inf*ones(4,1);
    b_U = [  -1 -1 -1 10]';
    c_L = zeros(3,1);
    c_U = zeros(3,1);
    x_opt = [ 0 1 1 1 0 1 1 0, 0.97, 0.9925, 0.98];
    f_opt = -0.94347;
    IntVars = 1:8;
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
elseif P == 418 % Problem from old glc_prob - NOT USED
    Name = 'Floudas-Pardalos 12.2 TP 1 IntV last';
    % Kocis and Grossmann (1988)
    x_L = [0    0  0 0 0 ]';
    x_U = [1.6 1.5 1 1 1 ]'; % Bnd x_1 from 3rd ineq., Bnd x_2 from 4th ineq
    A = [1 0 1 0 0;0 1.333 0 1 0; 0 0 -1 -1 1];
    b_L=[-inf -inf -inf]';
    b_U=[1.6 3 0];
    c_L = [1.25  3 ]';
    c_U = [1.25  3 ]';
    x_opt = [1.12 1.31 0 1 1];
    f_opt = 7.6672;
    IntVars = [3 4 5];
    x_min = x_L;
    x_max = x_U;
    x_0 = (x_L+x_U)/2;
    x_0(IntVars) = x_L(IntVars);
else
    error('glcIP_prob: Illegal problem number');
end

% Set x_0 to zeros
if isempty(x_0)
    x_0=zeros(length(x_L),1);
end

% Define the Prob
c  = 'glcIP_c';
if isempty(c_L) & isempty(c_U)
    c  = [];
end

Prob = glcAssign('glcIP_f', x_L, x_U, Name, A, b_L, b_U, ...
    c, c_L, c_U, x_0, IntVars, [], [], [], ...
    f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;
Prob.user = user;

% MODIFICATION LOG
%
% 080925  nhq  File created. IP-Problems merged from glc_prob and minlp_prob
% 080926  nhq  Added IntVars to problems 3,4,10 and 14 (formerly not IP's)
% 080929  nhq  File re-created. Names preserved from the old files.
% 080929  hkh  Revise comments
