% con_prob: Defines constrained nonlinear problems
%
% function [probList, Prob] = con_prob(P);
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

function [probList, Prob] = con_prob(P, varargin)

global Baux B1 B2 B3 B4 B5 B6 mass_vector gravity mass_matrix gravity_c

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Exponential-1. 2 inequalities'...
    ,'Exponential-2. 2 ineq. + bounds'...
    ,'Exponential-3. 2 ineq. + 1 eq.'...
    ,'Powell 1969'...
    ,'Fletcher 12.19'...
    ,'Circle-Triangle'...
    ,'Chvatal'...
    ,'Schittkowski 14'...
    ,'Schittkowski 24'...
    ,'Schittkowski 66'...
    ,'Fletcher page 279'...
    ,'ABB Robotics nonlinear problem'...
    ,'Hock-Shittkowski 375'...
    ,'DAS 2'...
    ,'Entropy problem'...
    ,'SOCS 6.3 Example 7.3.1'...
    ,'BMI rewritten as NLP'...
    ); % CHANGE: MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

HessPattern = []; ConsPattern = []; c_L = []; c_U = [];
pSepFunc = []; uP = []; f_Low = []; x_opt = []; f_opt = [];
A = []; x_min = []; x_max = [];

if P == 1
    % Exponential problem 1. 2 inequalities
    Name='Exponential-1. 2 inequalities';
    c_L=zeros(2,1);
    b_L=[]; b_U=[]; c_U=[];
    A=[];
    x_0=[-1;1];
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=0;
    x_min=[-10;0];
    x_max=[2;3];
    x_opt=[(-sqrt(449)-17)/4 (sqrt(449)-17)/4];
    f_opt=0.023550379624175;
elseif P == 2
    % Exponential problem 2. 2 inequalities + simple bounds 0 , x_L=0
    Name='Exponential-2. 2 ineq. + bounds';
    c_L=zeros(2,1);             %NHQ, 2nd constrint redundant.
    b_L=[]; b_U=[]; c_U=[];
    A=[];
    x_0=[-1;1];
    x_L=[0;0];
    x_U=[10;10];
    f_Low=0;
    x_min=[-2;-2];
    x_max=[3;3];
    x_opt=[0 1.5];
    f_opt=8.5;
elseif P == 3
    % Exponential problem 3. 2 inequalities + 1 linear equality. Not x_L >=0
    Name='Exponential-3. 2 ineq. + 1 eq.';
    % This problem probably have several minima.
    A=[1 1];
    b_L=0;
    b_U=0;
    c_L=zeros(2,1);
    c_U=[];
    x_0=[-5;5];
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=0;
    x_min=[-2;-2];
    x_max=[4;4];
    x_opt= [-3.162278, 3.162278]; % local min: [-1.224745, 1.224745];
    f_opt= 1.156627;              %  w. value: 1.8951
elseif P == 4
    Name='Powell 1969. Fletcher 12.10 page 327.';
    % 3 nonlinear equalities
    c_L=[10 0 -1]';
    c_U=[10 0 -1]';
    b_L=[]; b_U=[];
    A=[];
    x_0=[-2 2 2 -1 -1]';
    f_Low=0;
    x_L=-10*ones(5,1);
    x_U=10*ones(5,1);
    x_min=-3*ones(5,1);
    x_max=3*ones(5,1);
    x_opt=[-1.717108 1.595669 1.827311 -0.763647 -0.763647];
    f_opt=0.053949848364;
    ConsPattern=[ones(5,1),zeros(5,2)]';
    ConsPattern(2,2:5)=ones(1,4);
    ConsPattern(3,1:2)=ones(1,2);
elseif P == 5
    Name='Fletcher 12.19';
    % Page 329. One nonlinear inequality constraint
    c_L=zeros(1,1);
    b_L=[]; b_U=[]; c_U=[];
    A=[];
    x_0=[0;0];
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=-10;
    x_min=[-2;-2];
    x_max=[2;2];
    x_opt=[-0.5 0.25];
    f_opt=-0.25;
elseif P == 6
    Name='Circle-Triangle. Fletcher 12.21 page 330.';
    A=[1 0 -1  0
        0 1  0 -1
        0 0  1 -1
        ];
    b_L=[1 1 0]';
    b_U=inf*ones(3,1);
    c_L=-1;
    c_U=-1;
    x_0=[1 1 1 1]';
    % The other constraints imply x(i) >= 2 for i=1,2  and x(i)>=1 for i=3,4
    x_L=[-Inf -Inf -Inf 1]';
    x_U=100*ones(4,1);
    f_Low=0;
    x_min=[2 2 1 1]';
    x_max=[4 4 3 1.5];
    x_opt=[3.414213 3.414214 2.414213 1];
    f_opt=11.656854;
elseif P == 7
    Name='Chvatal';
    m = 5;
    uP(1)=m;
    A=[];
    b_L=[]; b_U=[]; c_U=[];
    c_L=zeros(m,1);
    x_0=zeros(m,1);
    x_L=zeros(m,1);
    x_U=10*ones(m,1);
    f_Low=-1E200;
    x_min=-1*ones(m,1);
    x_max=10*ones(m,1);
    x_opt=[];
    f_opt=[];
elseif P == 8
    % Schittkowski 14. Bracken, McCormick, Himmelblau. Start (2,2), f=1
    Name='Schittkowski 14';
    A=[1 -2];
    b_L=-1;
    b_U=-1;
    c_L=-1;
    c_U=[];
    x_0=[2;2];
    x_opt=[0.5*(sqrt(7)-1), 0.25*(sqrt(7)+1)];
    f_opt=9-2.875*sqrt(7);
%     x_opt=[0.5*(sqrt(13)-1), 0.25*(sqrt(13)+1)];
%     f_opt=(87-23*sqrt(13))/8;
    x_L=[-5;-5];
    x_U=[5;5];
    f_Low=0;
    x_min=[0;0];
    x_max=[2;2];
elseif P == 9
    % Schittkowski 24. Betts and Box. Start in (1,.5), f= -.01336459
    Name='Schittkowski 24';
    % 3 linear inequalities     %NHQ, 2nd constraint redundant
    r=sqrt(3);
    A=[1/r -1; 1 r; -1 -r];
    b_L=[0 0 -6]';
    b_U=inf*ones(3,1);
    c_L=[]; c_U=[];
    x_0=[1;0.5];
    x_L=[0;0];
    x_U=[5;5];
    f_Low=-10;
    x_min=[0;0];
    x_max=[4;2];
    x_opt=[3,sqrt(3)];
    f_opt=-1;
elseif P == 10
    % Schittkowski 66. Eckhardt. Start in (0,1.05,2.9). f=.58
    Name='Schittkowski 66. Eckhardt';
    % 2 nonlinear inequalities
    A=[];
    b_L=[]; b_U=[]; c_U=[];
    c_L=zeros(2,1);
    x_0=[0;1.05;2.9];
    x_L=[0;0;0];
    x_U=[100;100;10];
    f_Low=0;
    x_min=[0;0;0];
    x_max=[1;2;4];
    x_opt=[0.1841264879  1.202167873  3.327322322];
    f_opt=.5181632741;
    ConsPattern=[ 1 0; 1 1; 0 1]';
elseif P == 11
    %Fletcher page 279. Slow convergence for penalty problems
    Name='Fletcher page 279';
    A=[];
    b_L=[]; b_U=[];
    c_L=1;
    c_U=1;
    x_0=[0.0;0.0];
    x_L=[0;0];
    x_U=[2;2];
    f_Low=-10;
    x_min=[0;0];
    x_max=[2;2];
    x_opt=[1/sqrt(2),1/sqrt(2)];
    f_opt=-sqrt(2);
elseif P == 12
    % Nonlinear maximation for ABB Robotics, project spring 1996
    %    x(1:6) is the angle speed of each axis
    %    x(7:12) is the angle acceleration of each axis
    %    The speed is limited to -beta <= x(1:6) <= beta
    %
    % nonlinear constraints (quadratic):
    %    -alfa_i <= tau_i <= alfa_i
    alfa = [1650.77 1650.77 1650.77 495.6 544.3 190.2]';
    A=[];
    b_L=[]; b_U=[];
    c_L=-alfa;
    c_U=alfa;
    Name='Nonlinear maximization for ABB Robotics spring 1996';
    % Compute the problem variables using the following long code segment
    mass_vector=[3.1707408476445953 -0.32560151929832037  ...
        -0.32961783688790292 0.56908484945806381 2.51798528848921795E-2 ...
        0.12813021332566174];  % Row vector
    gravity=3.3434378370656517;
    gravity_c=[0 222.2222222222222 120.35527085274913
        4.9049995879799999 -11.787369547073027 3.5889336003556438]';
    mass_matrix=zeros(6,6);
    mass_matrix(1,:)=[128.9385764629213     -0.44499996261999997 ...
        5.80847885021473949E-3 2.6696882059615055  ...
        8.08084725502147411E-2 0.69481492454599469];
    mass_matrix(2,1:3)=[-0.44499996261999997 ...
        152.74598600000002  18.07411490421606];
    mass_matrix(2,5:6)=[-1.9993633957839383 0.30332922273072682];
    mass_matrix(3,:)=[5.80847885021473949E-3 18.07411490421606 ...
        83.640134496465876 -0.48980171872307404  ...
        1.6327180389337244 -0.36234174986602158];
    mass_matrix(4,1)=2.6696882059615055;
    mass_matrix(4,3:6)=[-0.48980171872307404 0.57350419891758142  ...
        -4.98017556830740726E-2 0.25365140250196505];
    mass_matrix(5,1)=0.42821593482321207;
    mass_matrix(5,2)=-1.8476987844185748;
    mass_matrix(5,3)=1.4515471640007136;
    mass_matrix(5,4)=7.70239455679084528E-2;
    mass_matrix(5,5)=0.70557835574547523;
    mass_matrix(5,6)=2.21284246518420463E-2;
    mass_matrix(6,1)=0.69481492454599469;
    mass_matrix(6,2)=0.30332922273072682;
    mass_matrix(6,3)=-0.36234174986602158;
    mass_matrix(6,4)=0.25365140250196505;
    mass_matrix(6,5)=-9.15215773121868525E-2;
    mass_matrix(6,6)=0.2273000039280578;
    Baux=zeros(6,6);
    Baux(1,1)=0.34002976897394821;
    Baux(1,2)=2.2644270560700002;
    Baux(1,3)=-0.12837583566546826;
    Baux(1,4)=5.91266598756878189E-2;
    Baux(1,5)=-0.51002084511546819;
    Baux(1,6)=8.08084725502147272E-2;
    Baux(2:6,1)=[2.2644270560700002; -0.12837583566546826;
        5.91266598756878189E-2; -0.51002084511546819; 8.08084725502147272E-2];
    Baux(2,2)=-0.30332922273072682;
    Baux(3,3:6)=[0.36234174986602158 -6.34008466993761371E-2 ...
        9.15215773121868525E-2 -8.66500084892087297E-2];
    Baux(4,3:6)=[-6.34008466993761371E-2 -3.64393902836788533E-2 ...
        0.35949593290521931 -4.98017556830740588E-2];
    Baux(5,3:6)=[9.15215773121868525E-2 0.35949593290521931  ...
        9.15215773121868525E-2 -8.66500084892087713E-2];
    Baux(6,3:6)=[-8.66500084892087297E-2 -4.98017556830740588E-2 ...
        -8.66500084892087713E-2 -9.15215773121868525E-2];
    B1=zeros(6,6);
    B1(1,2)=74.094318342615011;
    B1(1,3)=8.0758265886814726;
    B1(1,4)=8.08084725502147411E-2;
    B1(1,5)=-2.6696882059615055;
    B1(1,6)=0.42721155936915406;
    B1(2,1)=74.094318342615011;
    B1(3,1)=8.0758265886814726;
    B1(4,1)=8.08084725502147411E-2;
    B1(5,1)=-2.6696882059615055;
    B1(6,1)=0.42721155936915406;
    B1(3,3)=0.48980171872307404;
    B1(3,4)=-8.66500084892087574E-2;
    B1(3,5)=4.98017556830740726E-2;
    B1(3,6)=-6.34008466993760261E-2;
    B1(4,3)=-8.66500084892087574E-2;
    B1(4,4)=-0.58980171032307416;
    B1(4,5)=1.7597329922026472;
    B1(4,6)=-0.43155134628613506;
    B1(5,3)=4.98017556830740726E-2;
    B1(5,4)=1.7597329922026472;
    B1(5,5)=4.98017556830740726E-2;
    B1(5,6)=-6.3400846699376151E-2;
    B1(6,3)=-6.34008466993760261E-2;
    B1(6,4)=-0.43155134628613506;
    B1(6,5)=-6.3400846699376151E-2;
    B1(6,6)=-0.60696527402110201;
    B2=zeros(6,6);
    B2(1,1)=-74.094318342615011;
    B2(1,4)=-1.9993633957839383;
    B2(1,6)=-0.62299997241000027;
    B2(3,3)=-38.163918342615013;
    B2(3,4)=0.44499996261999997;
    B2(3,5)=-1.2321951426150108;
    B2(3,6)=0.32560151929832037;
    B2(4,1)=-1.9993633957839383;
    B2(4,3)=0.44499996261999997;
    B2(5,3)=-1.2321951426150108;
    B2(5,5)=-1.2321951426150108;
    B2(5,6)=0.32560151929832037;
    B2(6,1)=-0.62299997241000027;
    B2(6,3)=0.32560151929832037;
    B2(6,5)=0.32560151929832037;
    B2(6,6)=0.42466092868812394;
    B3=zeros(6,6);
    B3(1,1)=-8.0758265886814726;
    B3(1,4)=0.14988283158838311;
    B3(1,6)=8.52505604525888583E-2;
    B3(2,2)=38.163918342615013;
    B3(3,4)=-5.80847885021473949E-3;
    B3(3,5)=-2.1845719771709184;
    B3(3,6)=0.32961783688790292;
    B3(4,1)=0.14988283158838311;
    B3(4,3)=-5.80847885021473949E-3;
    B3(4,4)=-2.2203930608415194;
    B3(4,6)=-0.61387893508859048;
    B3(5,3)=-2.1845719771709184;
    B3(5,5)=-2.1845719771709184;
    B3(5,6)=0.32961783688790292;
    B3(6,1)=8.52505604525888583E-2;
    B3(6,3)=0.32961783688790292;
    B3(6,4)=-0.61387893508859048;
    B3(6,5)=0.32961783688790292;
    B3(6,6)=-0.50727846995863335;
    B4=zeros(6,6);
    B4(1,1)=-8.08084725502147411E-2;
    B4(1,2)=1.9993633957839383;
    B4(1,3)=-0.14988283158838311;
    B4(1,5)=-0.48685419042837269;
    B4(1,6)=5.50821870285079992E-2;
    B4(2,1)=1.9993633957839383;
    B4(2,2)=-0.44499996261999997;
    B4(3,1)=-0.14988283158838311;
    B4(3,3)=5.80847885021473949E-3;
    B4(3,5)=8.08084725502147411E-2;
    B4(3,6)=-5.90640043425956149E-2;
    B4(4,5)=0.24349442231358037;
    B4(4,6)=-3.3946806990795636E-2;
    B4(5,1)=-0.48685419042837269;
    B4(5,3)=8.08084725502147411E-2;
    B4(5,4)=0.24349442231358037;
    B4(5,5)=8.08084725502147411E-2;
    B4(5,6)=-5.90640043425956426E-2;
    B4(6,1)=5.50821870285079992E-2;
    B4(6,3)=-5.90640043425956149E-2;
    B4(6,4)=-3.3946806990795636E-2;
    B4(6,5)=-5.90640043425956426E-2;
    B4(6,6)=-6.23846544744587547E-2;
    B5=zeros(6,6);
    B5(1,1)=2.4560824262769283;
    B5(1,2)=0.31149998620500008;
    B5(1,3)=-4.26252802262944291E-2;
    B5(1,4)=0.4593130969141187;
    B5(1,5)=-9.51252779012944361E-2;
    B5(1,6)=0.19025055580258887;
    B5(2,1)=0.31149998620500008;
    B5(2,2)=1.0693943829658505;
    B5(3,1)=-4.26252802262944291E-2;
    B5(3,3)=2.019763058726967;
    B5(3,4)=-5.12764703789169163E-2;
    B5(3,5)=1.25899264424460897E-2;
    B5(3,6)=-2.5179852884892176E-2;
    B5(4,1)=0.4593130969141187;
    B5(4,3)=-5.12764703789169163E-2;
    B5(4,4)=-0.22652101881818254;
    B5(4,5)=-1.06051881570490975E-3;
    B5(4,6)=2.12103763140980561E-3;
    B5(5,1)=-9.51252779012944361E-2;
    B5(5,3)=1.25899264424460897E-2;
    B5(5,4)=-1.06051881570490975E-3;
    B5(5,5)=1.25899264424460897E-2;
    B5(5,6)=-2.51798528848921795E-2;
    B5(6,1)=0.19025055580258887;
    B5(6,3)=-2.5179852884892176E-2;
    B5(6,4)=2.12103763140980561E-3;
    B5(6,5)=-2.51798528848921795E-2;
    B5(6,6)=-0.12813021332566177;
    B6=zeros(6,6);
    B6(1,1)=-0.42721155936915406;
    B6(1,2)=0.62299997241000016;
    B6(1,3)=-8.52505604525888583E-2;
    B6(1,4)=-5.50821870285079992E-2;
    B6(1,5)=-0.19025055580258887;
    B6(2,1)=0.62299997241000016;
    B6(2,2)=-0.32560151929832037;
    B6(3,1)=-8.52505604525888583E-2;
    B6(3,3)=-0.32961783688790292;
    B6(3,4)=5.90640043425956426E-2;
    B6(3,5)=2.51798528848921795E-2;
    B6(3,6)=4.85722551762935077E-18;
    B6(4,1)=-5.50821870285079992E-2;
    B6(4,3)=5.90640043425956426E-2;
    B6(4,4)=3.3946806990795636E-2;
    B6(4,5)=-2.12103763140981949E-3;
    B6(5,1)=-0.19025055580258887;
    B6(5,3)=2.51798528848921795E-2;
    B6(5,4)=-2.12103763140981949E-3;
    B6(5,5)=2.51798528848921795E-2;
    B6(6,3)=4.85722551762935077E-18;

    % End of inirobot code segment

    beta=[2.6274; 2.0938; 2.0938; 3.9270; 4.2985; 5.7600];
    x_0=zeros(12,1);
    x_L=[-beta;-Inf*ones(6,1)];
    x_U=[beta; Inf*ones(6,1)];
    f_Low=-1000;
    x_min=x_L;
    x_max=x_U;
    x_opt=[];
    f_opt=[];
    HessPattern = [ones(6,6),zeros(6,6); zeros(6,12)];
elseif P == 13
    Name='Hock-Shittkowski 375';
    x_0   = ones(10,1);
    x_L   = [];
    x_U   = [];
    x_opt = [0.1064*ones(1,8), 2.843, -2.642];
    f_opt = -15.16;
    A    = ones(8,10);
    for i=1:8
        A(i,i) = 1/2;
    end
    b_L   = ones(8,1);
    b_U   = ones(8,1);
    c_L   = 4;
    c_U   = 4;
    x_min = -ones(10,1);
    x_max =  ones(10,1);
elseif P == 14
    Name='DAS2';
    % DAS ALOK EX. 2
    % Main1sub.for, row=1400
    x_0   = [1.4210526 0.9736842 0.1315789 1.5]';
    x_opt = [0.272727  2.090909  0 0.545455];
    f_opt = 4.18108504398827;
    x_L = [ 0   0   0   0 ]';
    x_U = [inf inf inf inf]';
    b_L=[-5;-4;1.5];
    b_U=[inf;inf;inf];
    A=[-1 -2 -1 -1;-3 -1 -2 1;0 1 4 0];
    c_L   = [];
    c_U   = [];
    x_min = -ones(4,1);
    x_max =  ones(4,1);
    pSepFunc=6;
elseif P == 15
    Name='Entropy problem';
    %-----------------------------------------------------------------------
    % 08 Oct 2002: Simple test program for pdco.m.
    %              "A" is an explicit sparse matrix, not a function.
    %              Michael Saunders, SOL, Stanford University.
    %-----------------------------------------------------------------------
    m = 50;
    uP(1)=m;
    n = 100;
    uP(2)=n;
    [A,b_U,x_L,x_U,d1,d2] = toydata( m,n );   % Private function below
    b_L = b_U;
    x_min = x_L;
    x_max = x_U;
    D  = sum(A,1);   D(D==0) = 1;
    D  = sparse( 1:n, 1:n, 1./D, n, n );
    A  = A*D;                                 % Normalize cols of A
    x_0 = ones(n,1)/n;
    Prob.SOL.d1 = d1;
    Prob.SOL.d2 = d2;
    options = pdcoSet;
    options.mu0       = 1e-1;  % 1.0 starts near central path
    options.LSQRatol1 = 1e-6;  % Let LSQR solve loosely to start with
    options.wait      = 0;     % Allow options to be reviewed before solve
    Prob.SOL.pdco     = options;
elseif P == 16
    % SOCS 6.3 Manual Example 7.3.1
    Name = 'SOCS 6.3 Example 7.3.1';
    x_L = [1e-4; 1e-4];
    x_U = [10; 10];
    b_L = [];
    b_U = [];
    c_L = 1;
    c_U = Inf;
    x_0 = [.5 ; 2];
    x_opt = [1 1];
    f_opt = 2;
elseif P==17
    % BMI problem bmi_prob(2) formulated as NLP.
    % c(1:2) expresses the bilinear inequalities,
    % and another linear semidefinite constraint
    % is contained within the simple bounds.
    Name = 'BMI rewritten as NLP';
    A = [1 1 1 0];
    b_L = -100;
    b_U = 100;
    x_L = [ 1   0  -inf -inf ]';
    x_U = [ inf inf inf inf ]';
    c_L = [0;0];
    c_U = [];
    x_0 = 2*ones(4,1);
else
    error('con_prob: Illegal problem number');
end

% Define the Prob
c  = 'con_c';
dc = 'con_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('con_f','con_g','con_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P  = P;
Prob.uP = uP;

function [A,b,bl,bu,d1,d2] = toydata( m,n )

% [A,b,bl,bu,d1,d2] = toydata( m,n );
% defines an m by n matrix A and rhs vector b,
% for use with pdco.m.
%
% 12 Feb 2001: First version of toydata.m.
% 30 Sep 2002: pdsco version modified for pdco.

rand('state',0);
density = 0.50;
rc      = 1e-1;
em      = ones(m,1);
en      = ones(n,1);
zn      = zeros(n,1);
if m==1
    A = sparse(ones(m,n));
else
    A = [sprand(m-1,n,density,rc)
        en'                     ];
end
x       = en/n;
b       = full(A*x);
bl      = zn;
bu      = en;
gamma   = 1e-4;
delta   = 1;          % Least squares
d1      = en*gamma;
d2      = em*delta;
d2(m)   = 1e-4;       % Make e'x = 1 satisfied more accurately.

% MODIFICATION LOG:
%
% 980826  hkh  Defining probType before call to ProbVarDef.
% 980922  hkh  Change name f_min to f_Low
% 981002  hkh  Change C_L to c_L for problem 5
% 981004  hkh  Add undefined b_L,b_U,c_L,c_U and A to all problems. Otherwise
%              problems if sending new values of these using structures
% 981006  hkh  Added call to checkuP
% 981011  hkh  Added call to checkuP
% 981011  hkh  Changed to use tomFiles for name definitions
% 981019  hkh  Update problem circle-triangle. Solution found.
% 981022  hkh  Set Name=[] if P is illegal
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 981106  mbk  Problem 13 added.
% 981116  hkh  Change f_Low for robot problem
% 011130  hkh  Change problem names, avoid too long names
% 030123  hkh  Add entropy problem
% 030126  hkh  New fields for pdco input, changes in problem 15
% 041102  ango Add SOCS 6.3 example
% 041117  med  xxx_prob removed and code added
% 050104  ango Added problem 17 - BMI as NLP
% 050406  hkh  Set bounds in comments for P 13,14, if using DIRECT
% 060116  med  x_L for 16 increased from 0 to 1e-4
% 080603  med  Switched to conAssign, cleaned
