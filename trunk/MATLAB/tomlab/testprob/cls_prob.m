% cls_prob: Defines constrained nonlinear least squares problems
%
% function [probList, Prob] = cls_prob(P);
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

function [probList, Prob] = cls_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'DIST'...
    ,'BAZA'...
    ,'GMW'...
    ,'PAW'...
    ,'DAS1'...
    ,'DAS2'...
    ,'DAS3'...
    ,'DAS4'...
    ,'Bob8'...
    ,'Bob9 - Infeasible'...
    ,'TP001'...
    ,'TP002'...
    ,'TP028'...
    ,'TP032'...
    ,'TP048'...
    ,'TP049'...
    ,'TP050'...
    ,'TP051'...
    ,'TP052'...
    ,'TP053'...
    ,'TP224'...
    ,'TP231'...
    ,'TP269'...
    ,'TP354'...
    ,'WrHo1'...
    ,'WrHo2'...
    ,'RELN'...
    ,'Constrained Walsh'...
    ,'EASY-TP14'...
    ,'EASY-TP6'...
    ,'EASY-TP43'...
    ,'EASY-TP57'...
    ,'EASY-TP327'...
    ,'EASY-TP355'...
    ,'EASY-GEO_PROB'...
    ,'EASY-PSS'...
    ,'EASY-4BAR_LNK'...
    ,'EASY-POL-APP'...
    ,'EASY-LIN_CMP1'...
    ,'EASY-DNS'...
    ,'EASY-RTD'...
    ,'EASY-RAT_APP'...
    ,'EASY-MIX_PAT1'...
    ,'EASY-MIX_PAT2'...
    ,'EASY-TREND'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

JacPattern = []; t = []; weightType = []; weightY = [];
SepAlg = []; f_Low = []; ConsPattern = []; f_opt = []; x_opt = [];
x_min = []; x_max = []; A = []; b_L = []; b_U = [];
c_L=[]; c_U=[];
uP = [];

if P==1
    Name='DIST';
    % A LEAST DISTANCE PROBLEM CONSTRUCTED BY PER L.
    % Main1sub.for, row=200
    m=2; n=2; % p=0; l=5;
    y=zeros(m,1);
    x_0=[10 6]';
    x_opt = [1.5 1.5];
    f_opt = 2.25;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
    A= [ 1.5 -1; 3 -1; 1 1; -0.5 1; -0.5 -1] ;
    b_L = [-2.5 -0.5 3 -1.5 -11]' ;
    b_U=Inf*ones(5,1);
elseif P==2
    Name='BAZA';
    % A LINEAR PROBLEM GIVEN BY BAAZARA AND SHETTY IN THEIR BOOK
    % SEE REPORT  "A WORKING ALGORITHM .......   PAGE NO. 43
    % Main1sub.for, row=300
    m=4; n=4; % p=0; l=4;
    y=zeros(m,1);
    x_0 = [-2 5 3 2 ]';
    x_opt = [0.952381 1 1.095238 0.809524];
    f_opt = 0.02380952380952381;
    x_L = [ 0   0   0   0 ]';
    x_U = [inf inf inf inf]';
    b_L=[-8;-2;-8;-6];
    b_U=Inf*ones(4,1);
    A=[-1 -1 -1 -1;-1 0 2 -4; -1 -1 0 0 ; 0 0 -1 -2];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==3
    Name='GMW';
    % A PROBLEM PROPOSED BY GILL,MURRAY,WRIGHT IN 'PRACTICAL
    % OPTIMIZATION' PAGE 171
    % Main1sub.for, row=1100
    m=2; n=2; % p=0; l=1;
    y=zeros(m,1);
    x_0 = [-3 2]';
    x_opt = [-2/3 -1/3];
    f_opt = 1/3;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=1;
    b_U=Inf;
    A=[-1 -1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==4
    Name='PAW';
    % A PROBLEM TO TEST WORKING SET STRATEGY
    % Main1sub.for, row=1200
    m=2; n=2; % p=0; l=2;
    y=zeros(m,1);
    x_0 = [0.4 0.5]';
    x_opt = [0.4 0.6];
    f_opt = 3.125;
    x_L = [ 0   0 ]';
    x_U = [inf inf]';
    b_L=[-1;-2.4];
    b_U=[inf;inf];
    A=[-1 -1;-3 -2];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==5
    Name='DAS1';
    % DAS ALOK EX. 1
    % Main1sub.for, row=1300
    m=2; % p=0; l=1;
    y=zeros(m,1);
    x_0 = [2 10]';
    x_opt = [2 0];
    f_opt = 0.02;
    x_L = [ 2 -50]';
    x_U = [50  50]';
    b_L=10;
    b_U=inf;
    A=[10 -1];
    x_min = [-1 -1]';
    x_max = [ 3 11]';
elseif P==6
    Name='DAS2';
    % DAS ALOK EX. 2
    % Main1sub.for, row=1400
    m=6; n=4; % p=0; l=7;
    y=zeros(m,1);
    x_0 = [1.4210526 0.9736842 0.1315789 1.5]';
    x_opt = [0.272727 2.090909 0 0.545455];
    f_opt = 4.18108504398827;
    x_L = [ 0   0   0   0 ]';
    x_U = [inf inf inf inf]';
    b_L=[-5;-4;1.5];
    b_U=[inf;inf;inf];
    A=[-1 -2 -1 -1;-3 -1 -2 1;0 1 4 0];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==7
    Name='DAS3';
    % DAS ALOK EX. 3
    % Main1sub.for, row=1500
    m=5; n=3; % p=0; l=4;
    y=zeros(m,1);
    x_0 = [0.5 0.5 0.5]';
    x_opt = [1.333333 0.777778 0.444444];
    f_opt = 77.20449172576832;
    x_L = [ 0   0   0 ]';
    x_U = [inf inf inf]';
    b_L=-3;
    b_U=inf;
    A=[-1 -1 -2];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==8
    Name='DAS4';
    % DAS ALOK EX. 4
    % Main1sub.for, row=1600
    m=11; n=10; l=8;
    y=zeros(m,1);
    x_0 = [0 0 15.0+1.0/3 -(15.0+1.0/3.0) 58 132 ...
        1.4285714 12.142856 -9.8181818 -30.545454]';
    x_opt = [2.629907 2.511192 10.000000 5.000000 2.203706 1.170634 ...
        1.494793 9.600977  8.790435 7.967652];
    f_opt = 103.3711489608575;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    A  = [ -4 -10   8  -3  -5  -5  -1   3
        -5   8  -2  -4  -8  -2  -2  -6
        0   0   0  -2  -1   0   0   0
        0   0   0   7   2   0   0   0
        0   0   0   0  -3  -14  0   0
        0   0   0   0   1   6   0   4
        3  17   0   0   0   0   0   0
        -9  -2   0   0   0   0   0   0
        0   0  -5   0   0   0   0 -12
        0   0   2   0   0   0   0   7 ]';
    b_L = -[105   0  12 138  46  42  20  96 ]';
    b_U =  inf*ones(l,1);
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==9
    Name='Bob8';
    % Main1sub.for, row=2500
    m=3; n=3; % p=0; l=3;
    y=zeros(m,1);
    x_0 = [0 0 0]';
    x_opt = [0.833333 0.666667 -0.166667];
    f_opt = 0.583333333332166;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[1.5;1;0.5];
    b_U=[inf;inf;inf];
    A=[1 1 0;1 0 -1;0 1 1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==10
    Name='Bob9 - Infeasible';
    % Main1sub.for, row=2600
    m=3; n=3; % p=0; l=3;
    y=zeros(m,1);
    x_0 = [0.81 0.61 -0.21]';
    x_opt = [];
    f_opt = [];
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[1.5;-1;-0.4];
    b_U=[inf;inf;inf];
    A=[1 1 0;-1 0 1;0 -1 -1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==11
    Name='TP001'; % Rosenbrock
    % Main2sub.for, row=100
    m=2; n=2; % p=0; l=1;
    y=zeros(m,1);
    x_0 = [-2 1]';
    x_opt = [1 1];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=-1.5;
    b_U=inf;
    A=[0 1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==12
    Name='TP002'; % Rosenbrock
    % Main2sub.for, row=200
    m=2; n=2; % p=0; l=1;
    y=zeros(m,1);
    x_0 = [-2 1]';
    x_opt = [1.224370749 1.5; -1.224371 1.5];
    f_opt = [0.025213093;2.4706146589945925];
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=1.5;
    b_U=inf;
    A=[0 1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==13
    Name='TP028';
    % Main2sub.for, row=1500
    m=2; n=3; % p=1; l=1;
    y=zeros(m,1);
    x_0 = [-4 1 1]';
    x_opt = [0.5 -0.5 0.5];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=1;
    b_U=1;
    A=[1 2 3];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==14
    Name='TP032';
    y=[0;0];
    x_0 = [0.1 0.7 0.2]';
    x_L = [0 0 0]';
    x_U = [100 100 100]';
    b_L=[1;3];
    b_U=[1;inf];
    A=[1 1 1;-1 6 4];
    x_min = [0 0 0]';
    x_max = [ 1  1  1]';
    x_opt = [0 0 1];
    f_opt = 0.5;
elseif P==15
    Name='TP048';
    % Main2sub.for, row=2100
    m=3; n=5; % p=2; l=2;
    y=zeros(m,1);
    x_0 = [3 5 -3 2 -2]';
    x_opt = [1 1 1 1 1];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[5;-3];
    b_U=[5;-3];
    A=[1 1 1 1 1;0 0 1 -2 -2];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==16
    Name='TP049';
    % Main2sub.for, row=2200
    m=4; n=5; % p=2; l=2;
    y=zeros(m,1);
    x_0 = [10 7 2 -3 0.8]';
    x_opt = [1 1 1 1 1];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[7;6];
    b_U=[7;6];
    A=[1 1 1 4 0;0 0 1 0 5];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==17
    Name='TP050';
    % Main2sub.for, row=2300
    m=4; n=5; % p=3; l=3;
    y=zeros(m,1);
    x_0 = [35 -31 11 5 -5]';
    x_opt = [1 1 1 1 1];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[6;6;6];
    b_U=[6;6;6];
    A=[1 2 3 0 0;0 1 2 3 0;0 0 1 2 3];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==18
    Name='TP051';
    % Main2sub.for, row=2400
    m=4; n=5; % p=3; l=3;
    y=zeros(m,1);
    x_0 = [2.5 0.5 2 -1 0.5]';
    x_opt = [1 1 1 1 1];
    f_opt = 0;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[4;0;0];
    b_U=[4;0;0];
    A=[1 3 0 0 0;0 0 1 1 -2;0 1 0 0 -1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==19
    Name='TP052';
    % Main2sub.for, row=2500
    m=4; n=5; % p=3; l=3;
    y=zeros(m,1);
    x_0 = [2 2 2 2 2]';
    x_opt = 1/349*[-33 11 180 -158 11];
    f_opt = 0.5*1859/349;
    x_L = -inf*ones(n,1);
    x_U =  inf*ones(n,1);
    b_L=[0;0;0];
    b_U=[0;0;0];
    A=[1 3 0 0 0;0 0 1 1 -2;0 1 0 0 -1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==20
    Name='TP053';
    % Main2sub.for, row=2600
    m=4; n=5; % p=3; l=13;
    y=zeros(m,1);
    x_0 = [2 2 2 2 2]';
    x_opt = 1/43*[-33 11 27 -5 11];
    f_opt = 0.5*176/43;
    x_L = -[10 10 10 10 10]';
    x_U =  [10 10 10 10 10]';
    b_L=[0;0;0];
    b_U=[0;0;0];
    A=[1 3 0 0 0;0 0 1 1 -2;0 1 0 0 -1];
    x_min = -ones(n,1);
    x_max =  ones(n,1);
elseif P==21
    Name='TP224';
    y=[0;0];
    x_0 = [0.1 0.1]';
    x_L = [0 0]';
    x_U = [6 6]';
    b_L=[0;-18;0;-8];
    b_U=inf*ones(4,1);
    A=[1 3;-1 -3;1 1;-1 -1];
    x_min = [-1 -1]';
    x_max = [ 7  7]';
    x_opt = [4 4];
    f_opt = 192;
elseif P==22
    Name='TP231'; % Rosenbrock
    y=[0;0];
    x_0 = [-1.2 1]';
    x_L = -[inf inf]';
    x_U =  [inf inf]';
    b_L=[-0.1;-0.1];
    b_U=[ inf; inf];
    A=[1/3 1;-1/3 1];
    x_min = [-1.5 -1]';
    x_max = [1.2 1.2]';
    x_opt = [1 1];
    f_opt = 0;
elseif P==23
    Name='TP269';
    y=[0;0;0;0];
    x_0 = [2 2 2 2 2]';
    x_L = -[inf inf inf inf inf]';
    x_U =  [inf inf inf inf inf]';
    b_L=[0;0;0];
    b_U=[0;0;0];
    A=[1 3 0 0 0;0 0 1 1 -2;0 1 0 0 -1];
    x_min = [0 0 0 0 0]';
    x_max = [1 1 1 1 1]';
    x_opt = [-0.7674 0.2558 0.6279 -0.1163 0.2558];
    f_opt = 0.5*4.09302;
elseif P==24
    Name='TP354';
    y=[0;0;0;0];
    x_0 = [3 -1  0  1]';
    x_L = -[inf inf inf inf]';
    x_U =  [20  20  20  20 ]';
    b_L=1;
    b_U=inf;
    A=[1 1 1 1];
    x_min = [0 0 0 0]';
    x_max = [1 1 1 1]';
    x_opt = [0.5034 -0.4557 0.2358 0.3064];
    f_opt = 0.5*0.113784;
elseif (P==25)|(P==26) % Osborne2, constraints from article of WRIGHT and HOLT
    if P==25
        x_0 = [1.3 0.65 0.65 0.7 0.6 3 5 7 2 4.5 5.5]'; % Original
        Name='WrHo1';
    else
        x_0 = [1.3 0.65 0.65 0.7 0.6 3 5 7 2 4.5 4.5]'; % Wright & Holt
        Name='WrHo2';
    end
    y = [1.366 1.191 1.112 1.013 0.991 0.885 0.831 0.847 0.786 0.725 0.746 0.679 ...
        0.608 0.655 0.615 0.606 0.602 0.626 0.651 0.724 0.649 0.649 0.694 0.644 ...
        0.624 0.661 0.612 0.558 0.533 0.495 0.500 0.423 0.395 0.375 0.372 0.391 ...
        0.396 0.405 0.428 0.429 0.523 0.562 0.607 0.653 0.672 0.708 0.633 0.668 ...
        0.645 0.632 0.591 0.559 0.597 0.625 0.739 0.710 0.729 0.720 0.636 0.581 ...
        0.428 0.292 0.162 0.098 0.054]';
    m=length(y);
    t=zeros(m,1);
    for i = 1 : m
        t(i) = (i-1)*0.1;
    end
    x_opt = [1.31 0.4315 0.6336 0.5993 0.7539 0.9056 1.3651 4.8248 2.3988 4.5689 5.6754];
    f_opt = 0.5*4.01683e-2;
    x_L = -inf*ones(length(x_0),1);
    x_U =  inf*ones(length(x_0),1);
    b_L=[6.270063;1.741584];
    b_U=[inf;inf];
    A=[1 2 3 4 0 0 0 0 0 0 0;1 0 1 0 0 0 0 0 0 0 0];
    x_min = -ones(length(x_0),1);
    x_max =  ones(length(x_0),1);
elseif P==27
    % Test of releasing more than one bound with variable dimension
    Name='RELN';
    % n problem variables, n >= 1 , default n = 10
    n     = 10; %Give problem dimension
    uP(1) = n;
    y    = zeros(n,1);
    x_0   = zeros(n,1);
    x_opt = 3.5*ones(1,n);
    f_opt = 0.5*n*0.25;
    x_L   = zeros(n,1);
    x_U   = 3.75*ones(n,1);
    b_L   = [-2*ones(n*(n-1),1);-3.5*n];
    b_U   = inf*ones(n*(n-1)+1,1);
    A     = zeros(n,n*(n-1)+1);
    con   = 0;
    for i = 1:n-1
        for j= i+1:n
            con=con+1;
            A(i,2*con-1) = -1;
            A(j,2*con-1) =  1;
            A(i,2*con)   =  1;
            A(j,2*con)   = -1;
        end
    end
    A(:,n*(n-1)+1) =  -ones(n,1);
    A=A';
    x_min =  -ones(n,1);
    x_max = 5*ones(n,1);
elseif P==28
    Name='Constrained Walsh';
    t=[ 2000; 5000; 10000; 20000; 30000; 50000];
    y=[0.9427; 0.8616; 0.7384; 0.5362; 0.3739; 0.3096];
    C = 96.05; %Give new C (Normal value 96.05)
    uP(1)=C;
    x_0=[-5;5];
    x_L=[-Inf; 1E-10];
    x_U=[1E-10; Inf];
    % Also constraint x(2) - t(i)*x(1) >= 0, to avoid log of negative
    A = [-50000 1];
    b_L=0;
    b_U=Inf;
    x_min=[-0.010 400]';
    x_max=[-0.008 650]';
elseif P==29
    Name='EASY-TP14';
    %Constrained least squares problem
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    m=2;
    y=zeros(m,1);
    x_0=[2 2]';
    x_opt = [0.8287566 0.91143783]'; % Estimated from EASY-FIT
    f_opt = 0.10330596e-03; % Estimated from EASY-FIT
    x_L = [0;0];
    x_U = [1000;1000];
    A= [1 -2];
    b_L =-1;
    b_U =[];
    c_L = -1;
    c_U = [];
elseif P==30
    Name='EASY-TP6';
    %Rosenbrock's banana function, Betts' formulation
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    m=1;
    y=zeros(m,1);
    x_0 = [100 100]';
    x_opt = [1.0000000 0.9999990]'; % estimated from EASY-FIT
    f_opt = 0.000000; % estimated from EASY-FIT
    x_L = [-1e4;-1e4];
    x_U = [1e5;1e5];
    x_min = [];
    x_max = [];
    c_L = [];
    c_U = 0;
elseif P==31
    Name='EASY-TP43';
    %Rosen-Suzuki test problem
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    m = 1;
    y=zeros(m, 1);
    x_0=zeros(4,1);
    x_opt = [-0.33537851e-9  1.000000  2.000000  -1.000000]'; % Estimated from EASY-FIT
    f_opt =0.4780000e3; % Estimated from eASY-FIT by KUHN-TUCKeR OPTIMALITY CRITeRION
    x_L = 1e3*(-ones(4,1));
    x_U = 1e5*(ones(4,1));
    x_min =[];
    x_max =[];
    c_L = [-8 -10 -5]';
    c_U = [];
elseif P==32
    Name='EASY-TP57';
    %Hock W., Schittkowski K.
    %Constrained exponential fit
    t = [8 8 10 10 10 10 12 12 12 12 14 14 14 16 16 16 18 18 20 20 20 22 22 22 24 24 24 26 26 26 28 28 30 30 30 32 32 34 36 36 38 38 40 42];
    y = [0.49 0.49 0.48 0.47 0.48 0.47 0.46 0.46 0.45 0.43 0.45 0.43 0.43 0.44 0.43 0.43 0.46 0.45 0.42 0.42 0.43 0.41 0.4 ...
        0.41 0.42 0.4 0.4 0.41 0.4 0.41 0.41 0.4 0.4 0.4 0.38 0.41 0.4 0.4 0.41 0.38 0.4 0.4 0.39 0.39];
    x_0 = [4.2000e-1 5.0000]';
    x_opt = [0.41995265 1.2848452]'; % Estimated from EASY-FIT
    f_opt = 0.14229835e-01; % Estimated from EASY-FIT
    x_L = [4.00e-1 -4.00e-1]';
    x_U = [1.00e2 2.56e1]';
    x_min =[];
    x_max =[];
    c_L = 0.09;
    c_U = [];
elseif P==33
    Name='EASY-TP327';
    t = [8 8 10 10 10 10 12 12 12 12 14 14 14 16 16 16 18 18 20 20 20 22 22 22 24 24 24 ...
        26 26 26 28 28 30 30 30 32 32 34 36 36 38 38 40 42];
    y = [0.49 0.49 0.48 0.47 0.48 0.47 0.43 0.46 0.46 0.45 0.45 0.43 0.43 0.44 0.43 0.43 0.46 ...
        0.45 0.42 0.42 0.43 0.4 0.41 0.41 0.4 0.4 0.42 0.41 0.4 0.41 0.41 0.4 0.4 0.4 0.38 0.41 0.4 ...
        0.4 0.41 0.38 0.4 0.4 0.39 0.39];
    x_0 = [4.2e-1 1.2848]';
    x_opt=[0.41995265 1.2848452]';
    f_opt = 0.14229835e-1;% Estimated from EASY-FIT
    x_L = [0 0]';
    x_U = [1.0e4 1.0e4]';
    x_min =[];
    x_max =[];
    c_L = 0.09;
    c_U = [];
elseif P==34
    Name='EASY-TP355';
    %Constrained least squares problem, four quadratic terms and local solutions
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    m=2;
    y=zeros(m,1);
    x_0=zeros(4,1);
    x_opt = [0.10000000 0.10000000 2.2000000 0.0000000]'; % Estimated from EASY-FIT
    f_opt = 0.61105000e2;  % Estimated from EASY-FIT by KUHN-TUCKeR OPTIMALITY CRITeRION
    x_L = zeros(4,1);
    x_U = 1e5*ones(4,1);
    x_min =[];
    x_max =[];
    A= [1 0 0 0;0 1 0 0];
    b_L = [0.1;0.1];
    b_U = [0.1;0.1];
    c_L =0;
    c_U =0;
elseif P==35
    Name='EASY-GEO_PROB';
    %Maximum distance from origin to intersection of ellipsoid with hyperboloid
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    m = 1;
    y=zeros(m,1);
    x_0 = [1.0 1.0 -1.0]';
    x_opt = [0.98842000 2.6736605 -1.8844639]'; % Estimated from EASY-FIT
    f_opt = 0.39005080e+04; % Estimated from EASY-FIT by KUHN-TUCKER OPTIMALITY CRITERION
    x_L = -1e3*(ones(3,1));
    x_U = 1e5*(ones(3,1));
    x_min = [];
    x_max = [];
    c_L = [14.92 2];
    c_U = [14.92 2]';
elseif P==36
    Name='EASY-PSS';
    %Primary and secondary stable model
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    t = [1 1.5 1.92 2.5 3.1];
    y = [33.5 52.4 58.6 72 74];
    x_0 = [1 1 1e1 1e1 -2.49e1]';
    x_opt = [1.8822582 23.188224 41.552647 82.999997 -24.900000]'; % Estimated from EASY-FIT
    f_opt = 0.30172637e-1; % Estimated from EASY-FIT
    x_L = [0 0 0 0 -2.49e1]';
    x_U = [1e3 1e3 1e3 1e3 -2.49e1]';
    x_min = [];
    x_max = [];
    A = [-1 0 1 0 0];
    b_L = -1e1;
    b_U = [];
elseif P==37
    Name='EASY-4BAR_LNK';
    %Design of a four bar linkage
    y = [1.00356431988 1.23508141645 1.48352986415 1.73197831184 1.96349540841 2.16230364633 2.31485457271 2.41075209405 ...
        2.44346095272 2.41075209408 2.31485457278 2.16230364642 1.96349540853 1.73197831197 1.48352986428 1.23508141657 ...
        1.00356431999 0.804756082062 0.652205155651 0.55630763429 0.523598775583 0.556307634196 0.652205155468 0.804756081847];
    for i=1:24
        t(i)=i;
    end
    x_0 = [2.0000e-1 1 3.0000e-1 4.0000e-1]';
    x_opt = [0.21420214 0.15404188 0.71695245e-1 -0.15001906e-1]'; % Estimated from EASY-FIT
    f_opt = 0.25839621; % Estimated from EASY-FIT
    x_L = [2.0000e-2 0 0 0]';
    x_U = [1 1.0000e6 2 1.0000e6]';
    x_min =[];
    x_max =[];
    A = [-1 0 0 1;-1 0 1 0];
    b_L = [1e-5;1e-5];
    b_U = [];
elseif P==38
    Name='EASY-POL-APP';
    %Hock W., Schittkowski K. (1981):
    %Polynomial approximation for computing axial forces
    t = [0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180];
    y = [0 0.0931 0.224 0.393713 0.579 0.7846 0.9535 1 1 1 1 1 0.95358 0.7846 0.579 0.3937 0.22425 0.0931 0];
    x_0 = [1 0 0 0 0 0 0 1 0 0 0 0 0 0]';
    x_opt = [-0.31615971e-03  0.72081494e-02  0.28957972e-03 -0.94846346e-05 ...
        0.34903272e-06 -0.55123896e-08  0.27420008e-10  0.10017580e1 ...
        -0.97884315e-2  0.14944132e-2 -0.69567717e-4  0.12024109e-5 ...
        -0.93163731e-8  0.27493302e-10]'; % Estimated from EASY-FIT
    f_opt =0.35361354e-3; % Estimated from EASY-FIT
    x_L = [-1.0000e3 -1.0000e4 -1.0000e4 -1.0000e4 -1.0000e4 -1.0000e4 -1.0000e4 ...
        -1.0000e3 -1.00e4 -1.0000e4 -1.00e4 -1.00e4 -1.00e4 -1.00e4]';
    x_U = 1e5*ones(14,1);
    x_min = [];
    x_max = [];
    c_L = zeros(19,1);
    c_U = [];
elseif  P==39
    Name='EASY-LIN_CMP1';
    %Heinzel G., Woloszczak R., Thomann P. (1993): TOPFIT 2.0:
    %Pharmaco-kinetic and Pharmakodynamic Data Analysis System, G. Fischer, Stuttgart
    %Linear compartments with bolus administration, single dose
    t = [0.25 0.5 1 1.5 2 4 6 12 24 72 120];
    y = [0 0 11.3931753437 15.5721363506 18.1009758214 26.0136307864 33.0546595141 50.4056725768 72.7820560242 97.5308438292 99.7760032059];
    x_0 = [50 -1 -1.5 1 1 1 0]';
    x_opt = [0.10000011e3 -0.29997517e1 -0.49999751e1  0.22747845e1 0.39997488  0.37499212  0.50002350]';
    f_opt = 0.18561345e-10;
    x_L = [0 -1.00e4 -1.00e4 1.00e-2 1.00e-2 1.00e-2 0]';
    x_U = [1.00e4 1.00e-3 1.00e-3 1.00e3 1.00e3 1.00e6 1.00e1]';
    x_min = [];
    x_max =[];
    A = [0 1 1 1 1 1 0;0 0 -1 1 0 0 0];
    b_L = [0;0.00001];
    b_U = [];
    c_L = 0;
    c_U = [];
elseif P==40
    Name='EASY-DNS';
    %Poeppe C., Pelliciari  C., Bachmann K. (1979):
    %Computer analysis of Feulgen hydrolysis kinetics,  Histochemistry, Vol. 60, 53-60
    %Feulgen-hydrolysis of DNS, biochemical reaction
    t = [6 12 18 24 30 36 42 48 54 60 66 72 78 84 90 96 102 108 114 120 126 132 138 144 150 156 162 168 174 180];
    y = [24.19 35.34  43.43 42.63 49.92 51.53 57.39 59.56 55.6 51.91 58.27 62.99 52.99 53.83 59.37 62.35 61.84 61.62 ...
        49.64 57.81 54.79 50.38 43.85 45.16 46.72 40.68 35.14 45.47 42.4 55.21];
    x_0=[8.77e2 9.12e-2 3.02e-3]';
    x_opt =[70.256884 0.50323123e-1 0.29789535e-2]';
    f_opt =0.49343685e-2;
    x_L = [0 0 0]';
    x_U = 1e5*ones(3,1);
    x_min = [];
    x_max = [];
    A = [0 1 -1];
    b_L = 0.001;
    b_U = [];
elseif P==41
    Name='EASY-RTD';
    %Schittkowski
    %Residence time distribution
    t = [6 12 18 24 30 36 42 48 54 60 66 72 78 84 90 96 102 108 114 120 126 132 138 144 150 156 162 168 174 180];
    y = [24.19 35.34  43.43 42.63 49.92 51.53 57.39 59.56 55.6 51.91 58.27 62.99 52.99 53.83 59.37 62.35 61.84 61.62 ...
        49.64 57.81 54.79 50.38 43.85 45.16 46.72 40.68 35.14 45.47 42.4 55.21];
    x_0 = [1e1 5]';
    x_opt = [59.630714 11.540072]';
    f_opt = 0.49343685e-2;
    x_L = [1e-10 1e-14]';
    x_U = [1e5 1e5]';
    x_min = [];
    x_max = [];
    A = [1 -1];
    b_L = 0.001;
    b_U = 0.001;
elseif P==42
    Name='EASY-RAT_APP';
    %Rational approximation with constraints
    %A general purpose algorithm for nonlinear least squares problems with nonlinear constraints,
    %Report UMINF-103.83,Institute of Information Processing, University of Umea, Umea, Sweden
    t = [0.0625 0.0714 0.0823 0.1 0.125 0.167 0.25 0.5 1 2 4];
    y = [0.0246 0.0235 0.0323 0.0342 0.0456 0.0627 0.0844 0.16 0.1735 0.1947 0.1957];
    x_0 = [2.50e-01 3.90e-1 4.15e-1 3.90e-1]';
    x_opt = [0.19226325  0.40401713  0.27497963  0.20678888]';
    f_opt = 0.20648571e-3;
    x_L = zeros(4,1);
    x_U = 1e5*ones(4,1);
    x_min =[];
    x_max =[];
    c_L =[0;0];
    c_U =[0.0246;0.1957];
elseif P==43
    Name='EASY-MIX_PAT1';
    %Schittkowski, Mixing pattern inside a polymerization reactor
    t = [0 0.004285 0.007829 0.009534 0.009813 0.011766 0.01348 0.011766 0.013673 0.011813 0.0122087 0.011464 0.01079 ...
        0.010673 0.011092 0.011441 0.010092 0.010232 0.00983 0.009929 0.00886 0.008092 0.008255 0.006883 0.007488 0.006534 ...
        0.005813 0.005302 0.004883 0.00393 0.003278 0.002302 0.002092 0.0021557 0.001418 0.000999 0.0006511 0.0003882];
    y = [0 0.003997 0.009097 0.012093 0.011737 0.012498 0.0108 0.0096 0.0097125  0.0097125 0.0088875 0.0082763 0.0073575 ...
        0.0066562 0.005265 0.0051 0.003937 0.003 0.0029437 0.00294 0.00195 0.0015375 0.0010337 0.00103 0.0008625 0.000525 0.00045 0.0002];
    x_0 = [3.0000e1 2.50000e1]';
    x_opt = [62.720945 6.0742287]';
    f_opt = 0.19176840e-2;
    x_L = [2.0000e1 1.0000e-3]';
    x_U = [1.0000e2 4.0000e1]';
    x_min = [];
    x_max = [];
    A = [-1 -1];
    b_L = -70;
    b_U = [];
elseif P==44
    Name='EASY-MIX_PAT2';
    %Schittkowski, Mixing pattern inside a polymerization reactor
    t = [0 5 7 9,9 11 13 16 19 22 24 26 28 30 33 36 39 42 45 48 53 58 68 83 98 113 128 143 158 173 188 218 248 278];
    y = [0 0.00784 0.015139 0.013934 0.012366 0.0115025 0.0133439 0.011025 0.010934 0.010859 0.011525 0.01122 0.010047 ...
        0.010525 0.009545 0.009331 0.00937 0.008561 0.007217 0.00695 0.005933 0.005728 0.005137 0.004478 0.002955 ...
        0.0024096 0.0021368 0.002068 0.0015912 0.001091 0.000568309 0.00034098 0.0001591];
    x_0 = [9.0000e1 3.0000e1 2.50000e1]';
    x_opt = [1.000000 64.905308 3.9995241]';
    f_opt = 0.60225436e-2;
    x_L = zeros(3,1);
    x_U = [1 1e5 1e5]';
    A = [0 -1 -1];
    b_L = -70;
    b_U = [];
elseif P==45
    Name='EASY-TREND';
    %Schittkowski Trend curve
    load('cls_trend.mat');
    x_0 = [8.358300 0.020450 0.000189 0.001046 0.598220 9.664700]';
    x_opt=[8.27 5.00e-3 -2.33e-4 -1.71e-4 8.00e-1 1.37e1]';
    f_opt = 0.19315285e-4;
    x_L = [0.0000 0.0000 -1.0000e2 -1.0000 1.0000e-7 -1.0000e2]';
    x_U = [1.0000e2 1.0000e2 1.0000e2 1.0000 1.0000e5 1.0000e5]';
    c_L = 0;
    c_U = [];
else
    error('cls_prob: Illegal problem number');
end

c = 'cls_c';
dc = 'cls_dc';

if isempty(c_L) & isempty(c_U)
    c = [];
    dc = [];
end

Prob = clsAssign('cls_r','cls_J', JacPattern, x_L, x_U, Name, x_0, ...
    y, t, weightType, weightY, SepAlg, f_Low, ...
    A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
    x_min, x_max, f_opt, x_opt);

if P ==45
    Prob.FUNCS.J = [];
end
Prob.P = P;
Prob.uP = uP;

% MODIFICATION LOG:
%
% 980826  hkh  Defining probType before call to ProbVarDef.
% 980922  hkh  Change name f_min to f_Low
% 981005  mbk  c_L=[]; c_U=[]; for all problems.
% 981006  hkh  Added call to checkuP
% 981011  hkh  Changed to use tomFiles for name definitions
% 981022  hkh  Set Name=[] if P is illegal
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 981210  hkh  Safeguard Walsh against division by 0 and log of negative numb
% 990623  hkh  Use clsVarDef and clsProbSet instead of ProbVarDef and ProbSet
% 000918  hkh  Safeguare Walsh bounding variables from zero
% 040517  med  Added 18 problems
% 041117  med  xxx_prob removed and code added
% 060814  med  FUNCS used for callbacks instead
% 080603  med  Switched to clsAssign, cleaned