% miqq_prob:
%
% Defines Mixed-Integer Quadratic Programming with Quadratic Constraints
% (MIQQ) problems
%
% function [probList, Prob] = miqq_prob(P);
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

function [probList, Prob] = miqq_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Fletcher EQP pg 231' ...
    ,'Fletcher IQP 10.4 pg 256' ...
    ,'Beales problem. Fletcher 8.25, pg 194' ...
    ,'Luenberger EQP pg 381'...
    ,'Fletcher 10.1 pg 255. EQP part.'...
    ,'Fletcher 10.1, page 255. IQP part.'...
    ,'Schittkowski IQP 21. Betts [8]'...
    ,'Fletcher EQP(n,K) 10.3'...
    ,'Pract.Opt. Ex.5.4 pg 171'...
    ,'Pract.Opt. Ex.5.5 pg 202'...
    ,'Bazaara EQP 11.14c. F neg.def.'...
    ,'Bazaara IQP 11.14c. F neg.def.'...;
    ,'Bazaara IQP 9.29b pg 405. F singular'...
    ,'Bunch and Kaufman Indefinite QP'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_opt = []; f_opt = []; x_min = []; x_max = [];
IntVars = []; x_L = []; x_U = [];

if P == 1
    % EQP-problem. Fletcher page 231.
    Name='Fletcher EQP pg 231';
    P1 = 3;
    x_0=[0 0 0]';
    x_L=[-10 -10 -10]';
    x_U=[10 10 10]';
    x_min=[0 0 -1]';
    x_max=[2 2 1]';
    c = zeros(3,1);
    if P1==1
        Name='Fletcher EQP pg 231';
        F = [2 0 0;0 2 0;0 0 2];
        A = [1 2 -1;1 -1 1];
        b_L = [4 -2]';
        b_U = b_L;
    elseif P1==2
        Name='Flet.EQP 1 Redundant';
        % Add 2*first row.
        F = [2 0 0;0 2 0;0 0 2];
        A = [1 2 -1;1 -1 1;2 4 -2];
        b_L = [4 -2 8 ]';
        b_U = b_L;
    elseif P1==3
        % Add 2*1st row, 4*first row and  2*2nd row as last
        Name='Flet.EQP 3 Redundant';
        F = [2 0 0;0 2 0;0 0 2];
        A = [1 2 -1;1 -1 1;2 4 -2;4 8 -4;1 -1 1];
        b_L = [4 -2 8 16 -2]';
        b_U = b_L;
    elseif P1==4
        % Change objective to f = x_1^2      + x_3^2
        Name='Flet.EQP F zero row';
        F = [2 0 0;0 0 0;0 0 2];
        A = [1 2 -1;1 -1 1];
        b_L = [4 -2]';
        b_U = b_L;
    elseif P1==5
        % Change objective to f = x_1^2+x_2^2+ 0.5E-17*x_3^2'...
        Name='Flet.EQP F(3,3) near 0';
        F = [2 0 0;0 2 0;0 0 1E-17]; % Matlabs QP fails for this
        A = [1 2 -1;1 -1 1];
        b_L = [4 -2]';
        b_U = b_L;
    end
    x_opt=[0.285714 1.428571 -0.857143];
    f_opt=2.8571428571428571428568;
elseif P==2
    % IQP-problem. Fletcher 10.4, page 256. Start in origo.
    % MATLAB OPT.TB FAILS in earlier version.
    %        gives (2,0), f=-2. Optimum (1.5,0.5), f=-2.75
    Name='Fletcher IQP 10.4 pg 256';
    beta = 1;
    F = beta*[2 -1;-1 2];
    c = [-3 0]';
    A = [-1 -1];
    b_L = -2;
    b_U = Inf;
    x_0=[0;0];
    x_min=[0;0];
    x_max=[5;5];
    x_L=[0;0];
    x_U=[100;100];
    x_opt=[1.5 0.5];
    f_opt=-2.7500000000000000;
elseif P==3
    % Beales LP problem. Fletcher 8.25, page 194.
    Name='Beales LP. Fletcher 8.25';
    F = zeros(4,4);
    c = [-0.75 20 -0.5 6]';
    A = [-0.25 8 1 -9;-0.5 12 0.5 -3];
    b_L = [0 0]';
    b_U = [Inf Inf]';
    x_0=zeros(4,1);	% Starting values for the optimization
    x_min=zeros(4,1);
    x_max=2*ones(4,1);
    x_L=[0;0;0;0];	% 4 used lower bounds for x
    x_U=[100;100;1;100];	% 1 used upper bound for x
    x_opt=[1 0 1 0];
    f_opt=-1.25;
elseif P==4
    Name='Luenberger EQP pg 381';
    n=10;
    v=zeros(n,1);
    for i = 1:n,v(i)=2*i;end
    F = diag(v);
    c = zeros(n,1);
    A = [1.5 1 1 0.5 0.5 0   0    0  0  0;
        0  0 0  0   0  2 -0.5 -0.5 1 -1;
        1  0 1  0   1  0   1    0  1  0;
        0  1 0  1   0  1   0    1  0  1];
    b_L = [5.5;2;10;15];
    b_U = b_L;
    x_0 = zeros(10,1);   % Not feasible
    x_L=-Inf*ones(10,1);
    x_U=Inf*ones(10,1);
    x_min=zeros(10,1);
    x_max=5*ones(10,1);
    x_opt=[-1.997878 2.664857 2.387961 3.622869 3.265128 2.865310 ...
        3.871821 3.158572 2.472968 2.688392];
    f_opt=502.43177928891447;
elseif P==5
    % Fletcher 10.1, page 255. EQP part of example. Same solution as IQP??
    Name='Fletcher EQP 10.1, pg 255';
    F = [ 3 -1  0
        -1  2 -1
        0 -1  1];
    c = ones(3,1);
    b_L = 4;
    b_U = b_L;
    A =  [ 1 2 1];
    x_0 = [0;0;0];
    x_L = [-Inf;-Inf;-Inf];
    x_U = [Inf;Inf;Inf];
    x_min = [0 0 0];
    x_max = [10 10 10];
    x_opt=[0.388889 1.222222 1.166667]';
    f_opt=3.277777777777777777;
elseif P==6
    % Fletcher 10.1, page 255.  IQP part of example. Same solution EQP??
    Name='IQP. Fletcher 10.1, pg 255';
    F = [ 3 -1  0
        -1  2 -1
        0 -1  1];
    c = ones(3,1);
    b_L = [ 4 0 0 0]';
    b_U = Inf*ones(4,1);
    A =  [ 1 2 1
        1 0 0
        0 1 0
        0 0 1 ];
    x_0 = [0;0;0];
    x_L = [-100;-100;-100];
    x_U = [100;100;100];
    x_min = [0 0 0];
    x_max = [10 10 10];
    x_opt=[0.388889 1.222222 1.166667];
    f_opt=3.27777777777777777;
elseif P==7
    % Schittkowski 21.QLR-T1-1. Betts[8].  IQP. Start with (-1,-1). Not feasible.
    % 1 constraint. 4 simple bounds. Optimum (2,0). f = 0.04.
    Name='IQP. Schittkowski 21';
    F = [ 0.02 0
        0    2];
    c = zeros(2,1);
    b_L =  10;
    b_U = Inf;
    A =  [ 10 -1];
    x_L=[2;-50];     % 2 used lower bounds for x
    x_U=[50;50];     % 2 used upper bounds for x
    x_0 = [ -1;-1]; % Start. Not feasible
    x_min = [-1 -1];
    x_max = [3 1];
    x_opt=[2 0];
    f_opt=0.04;
elseif P==8
    % Fletcher 10.3, page 255.  Min sum_i^n i x_i^2 s/t sum_i^n x_i = K
    Name='Fletcher EQP(n,K) 10.3';
    uP(2)=10;
    uP(1)=3;
    n=uP(1);
    F=diag(2*(1:n));
    c=zeros(n,1);
    A=ones(1,n);
    b_L=uP(2);
    b_U = b_L;
    x_0=zeros(n,1);
    x_L=-100*ones(n,1);
    x_U=100*ones(n,1);
    x_min=zeros(1,n);
    x_max=uP(2)*ones(1,n);
    if n==3 & b_L==10
        x_opt=[5.454545 2.727273 1.818182];
        f_opt=54.54545454545454;
    end
elseif P==9
    % Practical Optimization Ex.5.4, page 171
    % 1 inequality constraint.  Problem with Lagrange multipliers
    % Start at (-3,2). Optimum at (-2/3,-1/3);
    Name='Prac.Opt.Ex.5.4, pg 171. IQP';
    F = [ 2 0; 0 4];
    c = zeros(2,1);
    b_L = 1;
    b_U = Inf;
    A =  [ -1 -1];
    x_L=[-100;-100];
    x_U=[100;100];
    x_0=[-3;2];
    x_max=[0;3];
    x_min=[-4;-1];
    x_opt = [-2/3,-1/3];
    f_opt = 2/3;
elseif P==10
    % Practical Optimization Ex.5.5, page 202
    % Lower & Upper bounds.  Problem with 0 Lagrange multipliers
    % Start at (0,0). Optimum at (10,10);
    % IQP converges from all points. MATLAB:s QP fails.
    Name='Prac.Opt.Ex.5.5, pg 202. IQP. Zero mult';
    F = [ 0 -1; -1 0];
    c = zeros(2,1);
    A = [];
    b_L = [];
    b_U = [];
    x_L=[0;0];
    x_U=[10;10];
    x_0=[0;0];
    x_max=[4;4];
    x_min=[0;0];
    x_opt=[10 10];
    f_opt=-100;
elseif P==11
    % Bazaraa, Sherali, Shetty: NLP, Ex 11.14, page 541.
    % EQP + x>=0 simple bounds. F Negative definite
    % Start at origo (Not feasible). DIFFICULT!
    Name='Bazaraa EQP Ex.11.14c. F neg.def.';
    % Add simple bounds to A matrix
    F = [ -2 0 0 0 0; 0 -2 0 0 0;zeros(3,5)];
    c = [ 2;2;0;0;0];
    b_L = [ 4 4 8 0 0 0 0 0]';
    b_U = [ b_L(1:3);  Inf * ones(5,1)];
    A =  [-1  2 1 0 0; 1  1 0 1 0; 3 -2 0 0 1; eye(5)];
    x_0=[2;2;2;0;4]; %Infeas. qpSolve f=-2.88. qp f=-2.88. qpopt f=-2.88
    xopt1=[1.33333,2.66666,0,0,9.33333];
    xopt2=[3.2,0.8,5.6,0,0];
    xopt3=[0,2,0,2,12];
    x_opt=[xopt1;xopt2;xopt3];
    f_opt=[-0.888888;-2.88000000000000000;0];
elseif P==12
    % Bazaraa, Sherali, Shetty: NLP, Ex 11.14, page 541.
    % EQP + x>=0 simple bounds. F Negative definite. FORMULATED AS IQP, No slacks
    % Start at origo (Now feasible). DIFFICULT!
    % Local min at (0,0), f=0.
    Name='Bazaraa IQP 11.14c. F neg.def.';
    F = [ -2 0 ; 0 -2];
    c = [ 2;2];
    b_L = [ -4 -4 -8 ]';
    b_U = Inf*ones(3,1);
    A =  [1  -2 ; -1  -1 ; -3 2 ];
    % Start at (0,0) local min. (1,1) global max. (3.2,0.8) global min
    % Start with x_0 < (1,1) ==> convergence to local min at (0,0)
    x_L=[0;0];  % 2 used lower bounds
    x_U=[100;100];
    x_0=[0;0];
    x_min=[0 0];
    x_max=[5 5];
    x_opt = [3.2 0.8; 0 0];
    f_opt = [-2.88000000000;0];
elseif P==13
    % Bazaraa, Sherali, Shetty: IQP, Ex 9.29b, page 405.
    % Start at (0,3). Singular F
    % Local min at (3,0), f=15
    Name='Bazaraa Ex.9.29. IQP';
    F = [ 2 2 ; 2 2];
    c = [ 2;6];
    % Start at (0,0) local min. (1,1) global max. (3.2,0.8) global min
    x_L=[-1;-1];  % 2 unused lower bounds. Lower bounds implicit in A matrix
    A =  [1  1 ; 1  0 ; 0 1];
    b_L = [ 3 0 0 ]';
    b_U = Inf*ones(3,1);
    x_U=[100;100];
    x_opt=[3 0];
    f_opt=15;
    x_0=[0;3];
    x_min=[0 0];
    x_max=[5 5];
elseif P==14
    % Indefinite quadratic programming problem of Bunch and Kaufman,
    % A computational method for the indefinite quadratic programming
    % problem, Linear Algebra and its Applications, 34, 341-370 (1980).
    Name='Bunch and Kaufman Indefinite QP';
    x_min=[-2 -3 -4 -5 -6 -2 -3 -5];
    x_max=[2 3 2 1 0 7 8 9];
    c   = [  7   6   5   4   3   2   1   0 ]';
    A   = [ -1   1   0   0   0   0   0   0
        0  -1   1   0   0   0   0   0
        0   0  -1   1   0   0   0   0
        0   0   0  -1   1   0   0   0
        0   0   0   0  -1   1   0   0
        0   0   0   0   0  -1   1   0
        0   0   0   0   0   0  -1   1 ];
    x_L = [ -1 -2.1 -3.2 -4.3 -5.4 -6.5 -7.6 -8.7 ]';
    b_L = [ -1 -1.05 -1.1 -1.15 -1.2 -1.25 -1.3 ]';
    x_U = [  1   2    3    4    5    6    7    8]';
    b_U = [  inf inf inf inf inf inf inf ]';
    F   = [  1.69  1     2     3     4     5     6     7
        1     1.69  1     2     3     4     5     6
        2     1     1.69  1     2     3     4     5
        3     2     1     1.69  1     2     3     4
        4     3     2     1     1.69  1     2     3
        5     4     3     2     1     1.69  1     2
        6     5     4     3     2     1     1.69  1
        7     6     5     4     3     2     1     1.69 ];
    x_0 = [ -1    -2    -3    -4    -5    -6    -7    -8 ]';
    xopt1 = [ -1 -2 -3.05 -4.15 -5.3 6 7 8 ];
    xopt3 = [ 1 2 1.880144 .780144 -.369856 -1.569856 -2.819856 -4.119856 ];
    x_opt=  [xopt1;xopt3];
    f_opt=[-621.48782499999993;-131.77416786872979];
else
    error('miqq_prob: Illegal problem number');
end

Prob = miqqAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, [], ...
    IntVars, [], [], [], Name, [], [], x_min, x_max, f_opt, x_opt);
Prob.P = P;

% MODIFICATION LOG:
%
% 980826  hkh  Defining probType before call to ProbVarDef.
% 980922  hkh  Change name f_min to f_Low
% 980930  hkh  Problem 8. Conflict with uP that is defined for other problem.
%              When checking for nonempty, uP(1) is defined, and crash occurs.
%              nonempty.
% 981006  hkh  Added call to checkuP
% 981011  hkh  Changed to use tomFiles for name definitions
% 981022  hkh  Set Name=[] if P is illegal
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 011030  hkh  Changed problem names, avoiding too long names
% 041117  med  xxx_prob removed and code added
% 080603  med  Switched to conAssign, cleaned