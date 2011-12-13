% qp_prob:
%
% Defines Quadratic Programming problems
%
% function [probList, Prob] = qp_prob(P);
%
% INPUT:
%    P      Problem number
%           If isempty(P), return string matrix with problem names
%
% OUTPUT:
%    probList List of Problems
%    Prob   Problem Structure. See field list. Most fields are defined
%    here.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc. $Release: 6.2.0$
% Written June 1, 1999.   Last modified Jun 3, 2008.

function [probList, Prob] = qp_prob(P, varargin)

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
    ,'Bazaara IQP 11.14c. F neg.def.'...
    ,'Bazaara IQP 9.29b pg 405. F singular'...
    ,'Bunch and Kaufman Indefinite QP'...
    ,'cvxqp3_s'...
    ,'dpklo1'...
    ,'dualc5'...
    ,'gouldqp2'...
    ,'primalc5'...
    ,'qadlittl'...
    ,'qafiro'...
    ,'qbeaconf'...
    ,'qbore3d'...
    ,'qbrandy'...
    ,'qcapri'...
    ,'qe226'...
    ,'qgrow7'...
    ,'qisrael'...
    ,'qpcblend'...
    ,'qpcboei2'...
    ,'qpcstair'...
    ,'qrecipe'...
    ,'qsc205'...
    ,'qscagr25'...
    ,'qscfxm1'...
    ,'qscorpio'...
    ,'qscrs8'...
    ,'qscsd1'...
    ,'qsctap1'...
    ,'qshare1b'...
    ,'qstandat'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

f_Low = []; x_0 = []; x_min = []; x_max = []; x_opt = []; f_opt = [];

if P == 1
    % EQP-problem. Fletcher page 231.
    Name='Fletcher EQP pg 231';
    f_Low=-1E5;
    P1 = 1;
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
    x_0=zeros(4,1);
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
    for i = 1:n,v(i)=2*i; end
    F = diag(v);
    c = zeros(n,1);
    A = [1.5 1 1 0.5 0.5 0   0    0  0  0;
        0  0 0  0   0  2 -0.5 -0.5 1 -1;
        1  0 1  0   1  0   1    0  1  0;
        0  1 0  1   0  1   0    1  0  1];
    b_L = [5.5;2;10;15];
    b_U = b_L;
    x_0 = zeros(10,1);
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
    x_L=[2;-50]; % 2 used lower bounds for x
    x_U=[50;50]; % 2 used upper bounds for x
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
    x_L=[0;0];     % Lower bounds for x
    x_U=[10;10];   % Upper bounds for x
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
    % One Optimum at (1.33333,2.66666,0,0,9.33333) f=-0.888888; x_L=zeros(2,1)
    % Another Optimum at (3.2,0.8,5.6,0,0) f=-2.888888; x_L=zeros(3,1)
    % Another Optimum at (0,2,0,2,12) f=-0; x_L=zeros(5,1)
    % Add simple bounds to A matrix
    F = [ -2 0 0 0 0; 0 -2 0 0 0;zeros(3,5)];
    c = [ 2;2;0;0;0];
    b_L = [ 4 4 8 0 0 0 0 0]';
    b_U = [ b_L(1:3);  Inf * ones(5,1)];
    A =  [-1  2 1 0 0; 1  1 0 1 0; 3 -2 0 0 1; eye(5)];
    % Changing how many x_L and x_U set changes start solution
    % and also which optimum that is found
    % Explicitely setting lower bounds, even if they are present i A
    x_L=zeros(5,1);    % Lower bounds for x are set.
    x_U=100*ones(5,1); % Upper bounds for x
    %         x_0=zeros(5,1);  %Infeas. qpSolve f=-2.88. qp f=-2.88. qpopt f=-2.88
    %         x_0=[1;1;0;0;0]; %Infeas. qpSolve f=-2.88. qp f=-2.88. qpopt f=-2.88
    %         x_0=[1;1;3;2;7]; %Feas.   qpSolve f=0.     qp f=0.     qpopt f=-2.88
    %         x_0=[2;2;2;0;6]; %Feas.   qpSolve f=-0.88. qp f=0.     qpopt f=0
    %         x_0=[0;0;4;4;8]; %Feas.   qpSolve f=0.     qp f=0.     qpopt f=0
    %         x_0=[0;1;2;3;10];%Feas.   qpSolve f=0.     qp f=1.     qpopt f=0
    %         x_0=[.2;.1;4;3.7;7.6];%Feas.qpSolve f=0.   qp f=0.     qpopt f=0
    x_0=[2;2;2;0;4]; %Infeas. qpSolve f=-2.88. qp f=-2.88. qpopt f=-2.88
    x_min=zeros(5,1);
    x_max=5*ones(5,1);
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
    x_0 = [ -1  -2  -3  -4  -5  -6  -7  -8 ]';
    xopt1 = [ -1 -2 -3.05 -4.15 -5.3 6 7 8 ];
    xopt3 = [ 1 2 1.880144 .780144 -.369856 -1.569856 -2.819856 -4.119856 ];
    x_opt=  [xopt1;xopt3];
    f_opt=[-621.48782499999993;-131.77416786872979];
elseif P==15
    load qp_probmat P1;
    Name = P1.Name;
    F    = P1.F;
    c    = P1.c;
    A    = P1.A;
    b_L  = P1.b_L;
    b_U  = P1.b_U;
    x_L  = P1.x_L;
    x_U  = P1.x_U;
elseif P==16
    load qp_probmat P2;
    Name = P2.Name;
    F    = P2.F;
    c    = P2.c;
    A    = P2.A;
    b_L  = P2.b_L;
    b_U  = P2.b_U;
    x_L  = P2.x_L;
    x_U  = P2.x_U;
elseif P==17
    load qp_probmat P3;
    Name = P3.Name;
    F    = P3.F;
    c    = P3.c;
    A    = P3.A;
    b_L  = P3.b_L;
    b_U  = P3.b_U;
    x_L  = P3.x_L;
    x_U  = P3.x_U;
elseif P==18
    load qp_probmat P4;
    Name = P4.Name;
    F    = P4.F;
    c    = P4.c;
    A    = P4.A;
    b_L  = P4.b_L;
    b_U  = P4.b_U;
    x_L  = P4.x_L;
    x_U  = P4.x_U;
elseif P==19
    load qp_probmat P5;
    Name = P5.Name;
    F    = P5.F;
    c    = P5.c;
    A    = P5.A;
    b_L  = P5.b_L;
    b_U  = P5.b_U;
    x_L  = P5.x_L;
    x_U  = P5.x_U;
elseif P==20
    load qp_probmat P6;
    Name = P6.Name;
    F    = P6.F;
    c    = P6.c;
    A    = P6.A;
    b_L  = P6.b_L;
    b_U  = P6.b_U;
    x_L  = P6.x_L;
    x_U  = P6.x_U;
elseif P==21
    load qp_probmat P7;
    Name = P7.Name;
    F    = P7.F;
    c    = P7.c;
    A    = P7.A;
    b_L  = P7.b_L;
    b_U  = P7.b_U;
    x_L  = P7.x_L;
    x_U  = P7.x_U;
elseif P==22
    load qp_probmat P8;
    Name = P8.Name;
    F    = P8.F;
    c    = P8.c;
    A    = P8.A;
    b_L  = P8.b_L;
    b_U  = P8.b_U;
    x_L  = P8.x_L;
    x_U  = P8.x_U;
elseif P==23
    load qp_probmat P9;
    Name = P9.Name;
    F    = P9.F;
    c    = P9.c;
    A    = P9.A;
    b_L  = P9.b_L;
    b_U  = P9.b_U;
    x_L  = P9.x_L;
    x_U  = P9.x_U;
elseif P==24
    load qp_probmat P10;
    Name = P10.Name;
    F    = P10.F;
    c    = P10.c;
    A    = P10.A;
    b_L  = P10.b_L;
    b_U  = P10.b_U;
    x_L  = P10.x_L;
    x_U  = P10.x_U;
elseif P==25
    load qp_probmat P11;
    Name = P11.Name;
    F    = P11.F;
    c    = P11.c;
    A    = P11.A;
    b_L  = P11.b_L;
    b_U  = P11.b_U;
    x_L  = P11.x_L;
    x_U  = P11.x_U;
elseif P==26
    load qp_probmat P12;
    Name = P12.Name;
    F    = P12.F;
    c    = P12.c;
    A    = P12.A;
    b_L  = P12.b_L;
    b_U  = P12.b_U;
    x_L  = P12.x_L;
    x_U  = P12.x_U;
elseif P==27
    load qp_probmat P13;
    Name = P13.Name;
    F    = P13.F;
    c    = P13.c;
    A    = P13.A;
    b_L  = P13.b_L;
    b_U  = P13.b_U;
    x_L  = P13.x_L;
    x_U  = P13.x_U;
elseif P==28
    load qp_probmat P14;
    Name = P14.Name;
    F    = P14.F;
    c    = P14.c;
    A    = P14.A;
    b_L  = P14.b_L;
    b_U  = P14.b_U;
    x_L  = P14.x_L;
    x_U  = P14.x_U;
elseif P==29
    load qp_probmat P15;
    Name = P15.Name;
    F    = P15.F;
    c    = P15.c;
    A    = P15.A;
    b_L  = P15.b_L;
    b_U  = P15.b_U;
    x_L  = P15.x_L;
    x_U  = P15.x_U;
elseif P==30
    load qp_probmat P16;
    Name = P16.Name;
    F    = P16.F;
    c    = P16.c;
    A    = P16.A;
    b_L  = P16.b_L;
    b_U  = P16.b_U;
    x_L  = P16.x_L;
    x_U  = P16.x_U;
elseif P==31
    load qp_probmat P17;
    Name = P17.Name;
    F    = P17.F;
    c    = P17.c;
    A    = P17.A;
    b_L  = P17.b_L;
    b_U  = P17.b_U;
    x_L  = P17.x_L;
    x_U  = P17.x_U;
elseif P==32
    load qp_probmat P18;
    Name = P18.Name;
    F    = P18.F;
    c    = P18.c;
    A    = P18.A;
    b_L  = P18.b_L;
    b_U  = P18.b_U;
    x_L  = P18.x_L;
    x_U  = P18.x_U;
elseif P==33
    load qp_probmat P19;
    Name = P19.Name;
    F    = P19.F;
    c    = P19.c;
    A    = P19.A;
    b_L  = P19.b_L;
    b_U  = P19.b_U;
    x_L  = P19.x_L;
    x_U  = P19.x_U;
elseif P==34
    load qp_probmat P20;
    Name = P20.Name;
    F    = P20.F;
    c    = P20.c;
    A    = P20.A;
    b_L  = P20.b_L;
    b_U  = P20.b_U;
    x_L  = P20.x_L;
    x_U  = P20.x_U;
elseif P==35
    load qp_probmat P21;
    Name = P21.Name;
    F    = P21.F;
    c    = P21.c;
    A    = P21.A;
    b_L  = P21.b_L;
    b_U  = P21.b_U;
    x_L  = P21.x_L;
    x_U  = P21.x_U;
elseif P==36
    load qp_probmat P22;
    Name = P22.Name;
    F    = P22.F;
    c    = P22.c;
    A    = P22.A;
    b_L  = P22.b_L;
    b_U  = P22.b_U;
    x_L  = P22.x_L;
    x_U  = P22.x_U;
elseif P==37
    load qp_probmat P23;
    Name = P23.Name;
    F    = P23.F;
    c    = P23.c;
    A    = P23.A;
    b_L  = P23.b_L;
    b_U  = P23.b_U;
    x_L  = P23.x_L;
    x_U  = P23.x_U;
elseif P==38
    load qp_probmat P24;
    Name = P24.Name;
    F    = P24.F;
    c    = P24.c;
    A    = P24.A;
    b_L  = P24.b_L;
    b_U  = P24.b_U;
    x_L  = P24.x_L;
    x_U  = P24.x_U;
elseif P==39
    load qp_probmat P25;
    Name = P25.Name;
    F    = P25.F;
    c    = P25.c;
    A    = P25.A;
    b_L  = P25.b_L;
    b_U  = P25.b_U;
    x_L  = P25.x_L;
    x_U  = P25.x_U;
elseif P==40
    load qp_probmat P26;
    Name = P26.Name;
    F    = P26.F;
    c    = P26.c;
    A    = P26.A;
    b_L  = P26.b_L;
    b_U  = P26.b_U;
    x_L  = P26.x_L;
    x_U  = P26.x_U;
elseif P==41
    load qp_probmat P27;
    Name = P27.Name;
    F    = P27.F;
    c    = P27.c;
    A    = P27.A;
    b_L  = P27.b_L;
    b_U  = P27.b_U;
    x_L  = P27.x_L;
    x_U  = P27.x_U;
else
    error('qp_prob: Illegal problem number')
end

Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                  [], [], f_Low, x_min, x_max, f_opt, x_opt);
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
% 041115  med  Added 27 QP problems
% 041117  med  xxx_prob removed and code added
% 080603  med  Switched to conAssign, cleaned
