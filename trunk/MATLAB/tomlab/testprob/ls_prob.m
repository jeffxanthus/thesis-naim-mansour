% ls_prob: Defines nonlinear least squares problems
%
% function [probList, Prob] = ls_prob(P);
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

function [probList, Prob] = ls_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Powell '...
    ,'Walsh'...
    ,'Gisela'...
    ,'Circle Fitting'...
    ,'Population problem'...
    ,'Plasmid Stability n=2'...
    ,'Plasmid Stability n=3'...
    ,'Plasmid Stability n=3 (subst.)'...
    ,'Plasmid Stability n=3 (probability)'...
    ,'Parameterized test function (Huschens)'...
    ,'Signomial problem (rand)'...
    ,'Signomial problem (pseudorand)'...
    ,'Exponential problem (rand)'...
    ,'Exponential problem (pseudorand)'...
    ,'Trigonometric problem'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

JacPattern = []; t = []; weightType = []; weightY = [];
SepAlg = []; f_Low = []; f_opt = []; x_opt = []; uP = [];
b_L=[]; b_U=[]; A=[];
user = [];

if P==1
    Name='Powell';
    x_0=1; x_L=-50; x_U=50;
    x_min=-1; x_max=1;
    uP(1) = 0.1; %Give new LAMBDA
    y=[-1 1]';
elseif P==2
    Name='Walsh';
    t=[ 2000; 5000; 10000; 20000; 30000; 50000];
    y=[0.9427; 0.8616; 0.7384; 0.5362; 0.3739; 0.3096];
    C = 96.05; %Give new C (Normal value 96.05)
    uP(1)=C;
    x_0=[-0.01;400];
    x_L=[-Inf; 0];
    x_U=[   0; Inf];
    x_min=[-0.010 400]';
    x_max=[-0.008 650]';
elseif P==3
    Name='Gisela';
    t=[0.25; 0.5; 0.75; 1; 1.5; 2; 3; 4; 6; 8; 12; 24; 32; 48; 54; 72; 80;...
        96; 121; 144; 168; 192; 216; 246; 276; 324; 348; 386];
    y=[30.5; 44; 43; 41.5; 38.6; 38.6; 39; 41; 37; 37; 24; 32; 29; 23; 21;...
        19; 17; 14; 9.5; 8.5; 7; 6; 6; 4.5; 3.6; 3; 2.2; 1.6];
    K = 5; %Give new K (Normal value K=5)
    uP(1)=K;
    x_0=[6.8729,0.0108,0.1248]';
    x_L=[-Inf,-Inf,-Inf]';
    x_U=[Inf,Inf,Inf]';
    x_min=[6 0.005 0.1]';
    x_max=[7 0.05  0.2]';
elseif P==4
    [Name, x_0, x_L, x_U, f_Low, x_opt, f_opt, x_min, x_max, ...
        uP, t, y] = circleInit;
    y=y(:); % Store column by column
elseif P==5
    Name='Population problem';
    % Easy to change problem size.
    t=[-15; -10; -5; 0; 5; 10; 15];
    y=[60; 64; 71; 80; 90; 101; 116];
    x_0=[80,1]';
    x_L=[-Inf,-Inf]';
    x_U=[Inf,Inf]';
    x_min=[ 70 0 ]';
    x_max=[ 90 2 ]';
elseif P==6
    Name='Plasmid n=2';
    series = 1; %Choice of data series (1-6)
    uP(1)=series;
    [t, y, du_0, R_0, p_0, D] = plasinit(series);
    x_0 = [du_0 R_0]';
    x_L = 1E5*[-1 -1]';
    x_U = 1E5*[ 1  1]';
    x_min=[0 0]';
    x_max=[0.1 0.05]';
    uP(2) = p_0;
    uP(3) = D;
elseif P==7
    Name='Plasmid n=3';
    series = 1; %Choice of data series (1-6)
    uP(1)=series;
    [t, y, du_0, R_0, p_0, D] = plasinit(series);
    x_0 = [du_0 R_0 p_0]';
    x_L = 1e5*[-1 -1 0 ]';
    x_U = 1e5*[ 1  1 1 ]';
    x_min=[0    0    0  ]';
    x_max=[0.1  0.05 0.1]';
    uP(2) = p_0;
    uP(3) = D;
elseif P==8
    Name='Plasmid n=3 (subst.)';
    series = 1; %Choice of data series (1-6)
    uP(1)=series;
    [t, y, du_0, R_0, p_0, D] = plasinit(series);
    uminus_0 = D + du_0*(1-p_0);
    uplus_0  = D - du_0*p_0;
    x_0 = [uminus_0 uplus_0 R_0]';
    x_L =     [0 0 0]';
    x_U = 1E5*[1 1 1]';
    x_min=[0    0    0  ]';
    x_max=[0.1  0.05 0.1]';
    uP(2) = p_0;
    uP(3) = D;
elseif P==9
    Name='Plasmid n=3 (probability)';
    series = 1; %Choice of data series (1-6)
    uP(1)=series;
    [t, y, du_0, R_0, p_0, D] = plasinit(series);
    uminus_0 = D + du_0*(1-p_0);
    uplus_0  = D - du_0*p_0;
    x_0 = [uminus_0 uplus_0 R_0/uplus_0]';
    x_L = [ 0   0  0 ]';
    x_U = [1E5 1E5 1 ]';
    x_min=[0    0    0 ]';
    x_max=[0.1  0.05 1 ]';
    uP(2) = p_0;
    uP(3) = D;
elseif P==10
    Name='Parameterized test function (Huschens)';
    t=[];
    y=[2; 0; -1];
    phi = 0.5; %Give new phi
    uP(1)=phi;
    n = 2;
    x_0   = zeros(n,1);
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
elseif P==11
    Name='Signomial problem (rand)';
    % n variable, 1 <= n, default n=20
    % m variable, n <= m, default m=100
    % Standard test cases (n,m) = (10,50), (20,100)
    n = 20; %Give problem dimension
    uP(1)=n;
    m = 100; %Give m (m >= n)
    uP(2)=m;
    t     = [];
    y     = ones(m,1);
    AL = 0; AU = 3;
    PERC = min(90,100-200/n);
    l = 8;
    user.l = l;
    user.A = floor((AU-AL + 1E-2)*rand(m,n,l)+AL);
    user.A(100*rand(m,n,l) <= PERC) = 0;
    user.C = 10*(2*rand(m,l)-1);
    x_0   = zeros(n,1);
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
elseif P==12
    Name='Signomial problem (pseudorand)';
    % n variable, 1 <= n, default n=20
    % m variable, n <= m, default m=100
    % Standard test cases (n,m) = (10,50), (20,100)
    n = 20; %Give problem dimension
    uP(1)=n;
    m = 100; %Give m (m >= n)
    uP(2)=m;
    t     = [];
    y     = zeros(m,1);
    PERC = min(90,100-200/n);
    l = 8;
    user.l = l;
    user.A = zeros(m,n,l);
    user.C = zeros(m,l);
    na=1; nb=1000; nc=2000;
    CL = -100; CU = 100; AL = 0; AU = 3;
    % Compute the parameter values
    for i=1:m
        for k=1:l
            for j=1:n
                randv=pseudorand(na,nb,nc);
                na=randv(1);     nb=randv(2);     nc=randv(3);
                % Generate A values in [AL,AU]+1E-2
                user.A(i,j,k)=round((AU-AL+1.0e-2)*randv(4)+AL);
                randv=pseudorand(na,nb,nc);
                na=randv(1);     nb=randv(2);     nc=randv(3);
                % Set PERC of the A values as 0
                if 100*randv(4)<=PERC
                    user.A(i,j,k)=0;
                end
            end
            randv=pseudorand(na,nb,nc);
            na=randv(1);     nb=randv(2);     nc=randv(3);
            % Generate C values in [CL,CU]
            user.C(i,k)=(CU-CL)*randv(4)+CL;
        end
    end
    x_0=zeros(n,1);
    % Generate initial values in [-5,5]
    for i=1:n
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=10*randv(4)-5;
    end
    % Generate random y values in [-10,10]
    for i=1:m
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        y(i)=10*(randv(4)*2.0-1.0);
    end
    x_0=zeros(n,1);
    % Generate initial values in [-5,5]
    for i=1:n
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=10*randv(4)-5;
    end
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
elseif P==13
    Name='Exponential problem (rand)';
    % n variable, 1 <= n, default n=10
    % m variable, n <= m, default m=50
    % Standard test cases (n,m) = (10,50), (50,150)
    n = 10; %Give problem dimension
    uP(1)=n;
    m = 50; %Give m (m >= n)
    uP(2)=m;
    % Generate problem
    CL = -5; CU = 0; AL = -0.2; AU = 0.3;
    PERC = 50;
    l = 5;
    user.l = l;
    user.A = (AU-AL)*rand(m,n,l)+AL;
    user.A(100*rand(m,n,l) <= PERC) = 0;
    user.C = (CU-CL)*rand(m,l)+CL;
    t     = [];
    y     = ones(m,1);
    x_0   = zeros(n,1);
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
elseif P==14
    Name='Exponential problem (pseudorand)';
    % n variable, 1 <= n, default n=10
    % m variable, n <= m, default m=50
    % Standard test cases (n,m) = (10,50), (50,150)
    n = 10; %Give problem dimension
    uP(1)=n;
    m = 50; %Give m (m >= n)
    uP(2)=m;
    % Generate problem
    CL = -5; CU = 0; AL = -0.2; AU = 0.3; XL = -1; XU = 1;
    PERC = 50;
    l = 5;
    user.l = l;
    user.A = zeros(m,n,l);
    user.C = zeros(m,l);
    na=1; nb=1000; nc=2000;
    t     = [];
    % Compute the parameter values
    for i=1:m
        for k=1:l
            for j=1:n
                randv=pseudorand(na,nb,nc);
                na=randv(1);     nb=randv(2);     nc=randv(3);
                % Generate A values in [AL,AU]
                user.A(i,j,k)=(AU-AL)*randv(4)+AL;
                randv=pseudorand(na,nb,nc);
                na=randv(1);     nb=randv(2);     nc=randv(3);
                % Set PERC of the A values as 0
                if 100*randv(4)<=PERC
                    user.A(i,j,k)=0;
                end
            end
            randv=pseudorand(na,nb,nc);
            na=randv(1);     nb=randv(2);     nc=randv(3);
            % Generate C values in [CL,CU]
            user.C(i,k)=(CU-CL)*randv(4)+CL;
        end
    end
    x_0=zeros(n,1);
    % Generate initial values in [XL,XU] = [-1,1]
    for i=1:n
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=(XU-XL)*randv(4)+XL;
    end
    y=zeros(n,1);
    % Generate random y values in [-10,10]
    for i=1:m
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        y(i)=10*(randv(4)*2.-1.);
    end
    % Add the model values to the random y values in [-10,10]
    for i=1:m
        for k=1:l
            z=0;
            for j=1:n
                z=z+x_0(j)*user.A(i,j,k);
            end
            y(i)=y(i)+user.C(i,k)*exp(z);
        end
    end
    % Change the initial values in [XL,XU] = [-1,1]
    for i=1:n
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=x_0(i)+0.1*((XU-XL)*randv(4)+XL-x_0(i));
    end
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
elseif P==15
    Name='Trigonometric problem';
    % n variable, 1 <= n, default n=10
    % m variable, n <= m, default m=50
    % Standard test cases (n,m) = (10,50), (50,250)
    n = 10; %Give problem dimension
    uP(1)=n;
    m = 50; %Give m (m >= n)
    uP(2)=m;
    user.A = zeros(m,n);
    user.B = zeros(m,n);
    user.C = zeros(m,1);
    na=1; nb=1000; nc=2000;
    t     = [];
    y     = ones(m,1);
    x_0   = zeros(n,1);
    c1=zeros(n,1); s1=zeros(n,1);
    % Compute the parameter values
    for i=1:n
        for j=1:m
            randv=pseudorand(na,nb,nc);
            na=randv(1);     nb=randv(2);     nc=randv(3);
            user.A(j,i)=round(201.*randv(4)-100.);
            randv=pseudorand(na,nb,nc);
            na=randv(1);     nb=randv(2);     nc=randv(3);
            user.B(j,i)=round(201.*randv(4)-100.);
        end
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=pi*(2.*randv(4)-1.);
        c1(i)=cos(x_0(i));
        s1(i)=sin(x_0(i));
    end
    for i=1:m
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        % Alt 2
        % user.C(i)=10*(2.*randv(4)-1.);
        % Alt 1
        y(i)=10*(2.*randv(4)-1.);
        z = 0;
        for j=1:n
            % Alt 2
            % y(i)=y(i)+user.A(i,j)*s1(j)+user.B(i,j)*c1(j);
            % Alt 1
            z = z + user.A(i,j)*s1(j)+user.B(i,j)*c1(j);
        end
        user.C(i)=z;
    end
    for i=1:n
        randv=pseudorand(na,nb,nc);
        na=randv(1);     nb=randv(2);     nc=randv(3);
        x_0(i)=x_0(i)+pi*(randv(4)*2-1);
    end
    x_L   =-inf*ones(n,1);
    x_U   = inf*ones(n,1);
    x_min =-10*ones(n,1);
    x_max = 10*ones(n,1);
else
    error('ls_prob: Illegal problem number');
end

Prob = clsAssign('ls_r',[], JacPattern, x_L, x_U, Name, x_0, ...
    y, t, weightType, weightY, SepAlg, f_Low, ...
    A, b_L, b_U, [], [], [], [], [], ...
    x_min, x_max, f_opt, x_opt);

if any(P~=[11 12 13])
    Prob.FUNCS.J   = 'ls_J';
    Prob.FUNCS.d2r = 'ls_d2r';
end

Prob.P    = P;
Prob.uP   = uP;
Prob.user = user;

% =====================================================================
% plasinit.m
% =====================================================================
%
% function [t, y, du_0, R_0, p_0, D, nopuilr] = plasinit(series)
%
% plasinit is called by ls_prob to initiate parameters and
% determine starting values for du, R and p_ in the plasmid
% stability problem.
%
% INPUT PARAMETERS
% series    Measurement series for which starting values is to
%           be determined.
%
% OUTPUT PARAMETERS
% t         Time vector (in generations)
% y        Measured values of p_
% du_0      Growth rate difference (per generation)
% R_0       Rate of conversion from plasmid-bearing to plasmid-
%           free cells (segregational instability), (per generation)
% p_0       Value of p_ at t=0
% D         Dillution rate (per generation)
% nopuilr   Optimal number of points used in linear regression
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written Mar 17, 1998.  Last modified June 21, 1999.
%

function [t, y, du_0, R_0, p_0, D] = plasinit(series)

if nargin < 1
    series = 1;
end

switch series
    case 1 % Data980209.xls, D=0.22 1/h
        th = [0 16 28 50 76 100 121 158]';
        %tg = [0 5 8.9 15.8 24 31.7 38.4 50]';
        pplus  = 1e-2*[100 97 69 57 30 25 15  7]';
        pminus = 1e-2*[  0  3 31 43 70 75 85 93]';
        Dprime = 0.22;
    case 2 % Data980209.xls, D=0.16 1/h
        th = [0 19 39 64 88 115]';
        %tg = [0 4.4 9 14.8 20.3 26.5]';
        pplus  = 1e-2*[100 99 70 56 33 22]';
        pminus = 1e-2*[  0  1 30 44 67 78]';
        Dprime = 0.16;
    case 3 % Data980209.xls, D=0.38 1/h
        th = [0 24 50 72 100 149 172]';
        %tg = [0 13.15 27.4 39 55 81.7 94.3]';
        pplus  = 1e-2*[100 98 81 59 48 36 26]';
        pminus = 1e-2*[  0  2 19 41 52 64 74]';
        Dprime = 0.38;
    case 4 % Data980317.xls, D=0.33 1/h
        th = [3 19 30 43 52 71 76 93 101 116]';
        %tg = [1.5 9.6 14.4 20.6 25 34.1 36.5 44.7 48.5 55.7]';
        pplus  = 0.5*1e-2*[200 194 177 166 108 56 43 22 4 0]';
        pminus = 0.5*1e-2*[0 6 23 34 92 144 157 178 196 200]';
        Dprime = 0.33;
    case 5 % Data980317.xls, D=0.39 1/h
        th = [23 49 61 71 87 110 119 132 143 156]';
        %tg = [12.8 27.3 34 39.5 48.5 61.2 66.2 73.5 79.6 86.8]';
        pplus  = 0.5*1e-2*[200 119 198 179 151 56 34 22 12 12]';
        pminus = 0.5*1e-2*[0 3 2 21 49 144 166 178 74 188]';
        Dprime = 0.39;
    case 6 % Data980317.xls, D=0.25 1/h
        th = [50 71 96 117 140 165 191 207 218 235 258 267 280 291 304.5]';
        %tg = [15 21.5 29 35.4 42.5 50 58 62.7 66 71.1 78.1 81 85 88 92.2]';
        pplus  = 0.5*1e-2*[200 182 87 91 79 69 69 100 76 74 42 29 10 8 0]';
        pminus = 0.5*1e-2*[0 18 113 109 121 131 131 100 124 126 158 171 190 192 200]';
        Dprime = 0.25;
end
m   = length(th);
Td  = log(2)/Dprime;
t   = th./Td;
y  = pminus;
p_0 = pminus(1);
D   = Dprime/log(2);

% STARTING VALUES
% Case (a)
j = 0;
for i = 1:m % Determine pminus_a, the nonzero elements in pminus
    if pminus(i) > 1E-12
        j = j + 1;
        pminus_a(j) = pminus(i);
        t_a(j) = t(i);
    end
end
m_a = length(t_a);
f0_a=inf;
for i = 1 : m_a-1
    MC   = polyfit(t_a(i:m_a),log(pminus_a(i:m_a)),1);
    du = MC(1);
    R  = MC(1)*(exp(MC(2))-p_0);
    r = zeros(1,m);
    for j=1:m
        denom = ( p_0*du + R )*exp( (du+R)*t(j) ) - R*(1-p_0);
        numer = ( p_0*du + R )*exp( (du+R)*t(j) ) + du*(1-p_0);
        r(j) = denom/numer - y(j);
    end
    f0 = 0.5*r*r';
    if f0 < f0_a
        f0_a  = f0;
        du_a  = du;
        R_a   = R;
    end
end

% Case (b)
j = 0;
for i = 1:m % Determine pplus_b, the nonzero elements in pplus
    if pplus(i) > 1E-12
        j = j + 1;
        pplus_b(j) = pplus(i);
        t_b(j) = t(i);
    end
end
m_b = length(t_b);
f0_b=inf;
for i = 1 : m_b-1
    MC =  polyfit(t_b(i:m_b),log(pplus_b(i:m_b)),1);
    du =  MC(1)*(exp(-MC(2))-1);
    R  = -MC(1)*exp(-MC(2));
    r = zeros(1,m);
    for j=1:m
        denom = ( p_0*du + R )*exp( (du+R)*t(j) ) - R*(1-p_0);
        numer = ( p_0*du + R )*exp( (du+R)*t(j) ) + du*(1-p_0);
        r(j) = denom/numer - y(j);
    end
    f0 = 0.5*r*r';
    if f0 < f0_b
        f0_b  = f0;
        du_b  = du;
        R_b   = R;
    end
end

% Case (c)
j = 0;
for i = 1:m % Determine pminus_c, the elements in pminus < 1
    if pminus(i) < (1 - 1E-12)
        j = j + 1;
        pminus_c(j) = pminus(i);
        t_c(j) = t(i);
    end
end
m_c = length(t_c);
pstar = pminus_c(m_c);
f0_c=inf;
for i = 1:m_c-2
    MC =  polyfit(t_c(i:m_c-1),log(abs(pminus_c(i:m_c-1)-pstar)),1);
    du =  MC(1)/(1-pstar);
    R  = -MC(1)*pstar/(1-pstar);
    r = zeros(1,m);
    for j=1:m
        denom = ( p_0*du + R )*exp( (du+R)*t(j) ) - R*(1-p_0);
        numer = ( p_0*du + R )*exp( (du+R)*t(j) ) + du*(1-p_0);
        r(j) = denom/numer - y(j);
    end
    f0 = 0.5*r*r';
    if f0 < f0_c
        f0_c  = f0;
        du_c  = du;
        R_c   = R;
    end
end

% Determine best starting values, case (a), (b) or (c)
f0_abc  = [f0_a  f0_b  f0_c];
du_abc  = [du_a  du_b  du_c];
R_abc   = [R_a   R_b   R_c];

[f0_min idx] = min(f0_abc);
du_0    = du_abc(idx);
R_0     = R_abc(idx);

% =====================================================================
% circleinit.m
% =====================================================================
%  		circleInit.m
%
% function [Name, x_0, x_L, x_U, f_Low, x_opt, f_opt, x_min, x_max, ...
%           uP,  t, y] = circleInit(P, ask, uP);
%
% Defines a normally distributed point set around a circle
% with center c and radius r
%
% INPUT:  c       Center of circle
%         r       Radius of circle
%         seedval Initial random seed value(123456)
%         ask     If true ask for: # of points (100),
%                 std deviation of radius (1),
%
% OUTPUT: points  Points in the plane - matrix size(# of points, 2)
%         seedval Initial random seed value used
%
% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Written May 14, 1995. Last modified June 21, 1999.

function [Name, x_0, x_L, x_U, f_Low, x_opt, f_opt, x_min, x_max, ...
    uP,  t, y] = circleInit()

Name='Circle Fitting';
CentCi=[10;10];
radius=2;
seed_val=123456;

m=100; %Give number of random points
sigma=1; %Give standard deviation for radius

uP(1:6)=[CentCi(:);radius;seed_val;m;sigma];

rand('seed',seed_val);

theta=rand(m,1)*2*pi;
rdev=randn(m,1)*sigma;

y=zeros(m,2);
for i=1:m
    y(i,1)=(radius+rdev(i))*cos(theta(i))+CentCi(1);
    y(i,2)=(radius+rdev(i))*sin(theta(i))+CentCi(2);
end

x_L=[0,min(y)]';
x_U=[Inf,max(y)]';

% Estimate starting values for x:  x_0 = [ radius centre_x centre_y]
x_0 = sum(y)/m; % Mean point value is starting point for centre of circle
r = 0;
maxr=0;
for i = 1:m
    d = norm(y(i,:)-x_0);
    maxr = max(d,maxr);
    r = r + d;
end
x_0 = [r/m; x_0'];           % Mean radius is starting value
x_U(1)=maxr;                 % The worst radius possible
x_opt=[radius;CentCi(:)];    % Optimal point (in theory)
x_min=x_L;
x_max=x_U;
t=[];
f_Low=0;
f_opt=[];

function y=pseudorand(na,nb,nc)
% =====================================================================
% Initialize by setting na, nb, nc, in ranges (1,999),
% (1000,1999), (2000, 3136). Result in (0,1).
for i=1:4
    temp=nb+nc-3137;
    if temp<=0
        temp=temp+3137;
    end
    y(3)=nb;
    y(2)=na;
    y(1)=temp;
end
y(4)=temp/3137;

% MODIFICATION LOG
%
% 981006  hkh  Added call to checkuP. Modified for new TOMLAB
%              Rename from circinit to circleInit
% 981129  hkh  Missing ; on Itext defintion
% 030323  hkh  Added four problems
% 030526  ango Problem number P was accidentally printed, removed
% 041117  med  xxx_prob removed and code added
% 050302  hkh  Added problems 11-15, 12,14,15 similar to Wang,Li,Qi paper
% 050302  hkh  New sub function pseudorand used for problems 12,14,15
% 080603  med  Switched to clsAssign, cleaned