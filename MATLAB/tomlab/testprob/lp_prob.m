% lp_prob: Defines Linear Programming problems
%
% function [probList, Prob] = lp_prob(P);
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

function [probList, Prob] = lp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    '1st simple LP. OR course'...
    ,'Luenberger 2.6-2'...
    ,'Example LP-3 OA 5p'...
    ,'Luenberger 2.6-3'...
    ,'Luenberger 2.6-5'...
    ,'Marshall-Suurballe cycling example'...
    ,'1st simple LP. OR course (ineq.form)'...
    ,'Kuhn example.'...
    ,'Beale cycling example.'...
    ,'Winston Ch.6 Review 17. Dual feasible'...
    ,'Fletcher 8.21.3.9. Redundancy in constr.'...
    ,'Winston Ex. 4.12 B4. Max | |. Rewritten'...
    ,'adlittle'...
    ,'afiro'...
    ,'bandm'...
    ,'boeing2'...
    ,'bore3d'...
    ,'brandy'...
    ,'capri'...
    ,'e226'...
    ,'etamacro'...
    ,'farm'...
    ,'gams10a'...
    ,'gams30a'...
    ,'grow7'...
    ,'israel'...
    ,'kb2'...
    ,'kleemin8'...
    ,'lotfi'...
    ,'multi'...
    ,'nug05'...
    ,'orswq2'...
    ,'p0033'...
    ,'p0040'...
    ,'p0201'...
    ,'p0282'...
    ,'p0291'...
    ,'p0548'...
    ,'recipe'...
    ,'refine'...
    ,'sc105'...
    ,'sc205'...
    ,'sc50a'...
    ,'sc50b'...
    ,'scagr25'...
    ,'scagr7'...
    ,'scfxm1'...
    ,'scorpion'...
    ,'sctap1'...
    ,'share1b'...
    ,'share2b'...
    ,'stocfor1'...
    ,'vtp-base'...
    ,'zed'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

x_0 = []; f_opt = []; x_opt = []; f_Low = [];
x_min = []; x_max = [];

if P == 1
    Name='1st simple LP. OR course';
    % Problem formulated on standard form
    A   = [1 2 1 0
        4 1 0 1 ];
    b_U = [6;12];
    c   = [-7  -5  0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_opt=[18/7 12/7 0 0];
    x_U=[18/7 12/7 6 12]';
    f_opt=-26-4/7;
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=3*ones(n,1);
elseif P == 2
    Name='Luenberger 2.6-2';
    % Problem formulated on standard form
    % Shows that it is unnessary to add sum(x(i)=1,
    % because summing the two first constraints gives sum(x(i))=1
    A=[0.1 0.25 0.5 0.75 0.95; 0.9 0.75 0.5 0.25 0.05;ones(1,5) ];
    b_U=[0.3 0.7 1]';
    c   = [5  4 3 2 1.5]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    %B=~[0 1 1 0 0]';
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 3
    Name='Example LP-3 from OA 5p ';
    % Problem formulated on standard form
    A1=[ones(1,3),zeros(1,3);zeros(1,3),ones(1,3)];
    A2=[eye(3),eye(3);[6 0.4 0.26],zeros(1,3);zeros(1,3),[12/5 4/25 13/125]];
    A3=[12 7.5 10 0 0 0;0 0 0 24/5 3 4];
    A  = [A1,zeros(2,7);A2,eye(5),zeros(5,2);A3,zeros(2,5),-eye(2)];
    b_U= [10 25 12 15 25 10 8 90 95]';
    c  = [1500 2400 3000 1500 2400 3000,zeros(1,7)]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 4
    Name='Luenberger 2.6-3';
    % Problem formulated on standard form
    A  = [0.3 0 0 0.3 0 0;0 0.2 0 0 0.4 0; 0 0 0.3 0 0 0.2];
    b_U= [900000 800000 500000 ]';
    c  = [35 35 35 30 30 30]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 5
    Name='Luenberger 2.6-5';
    % Problem formulated on standard form, by adding two vars
    A   = [2 -2 -2 1; 1 -1 0 -1];
    b_U = [4 1]';
    c   = [1 -1 4 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 6
    Name='Marshall-Suurballe cycling example.';
    A   = -[-.5 5.5 2.5 -9  -1 0 0
        -.5 1.5 0.5 -1  0 -1 0
        -1 0   0    0  0 0 -1];
    b_U = [0 0 1]';
    c   = -[10 -57 -9 -24 0 0 0]';
    %B=~[0 0 0 0 1 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 7
    %Name='First LP. Operations Research 5p. (Inequality form)';
    Name='1st simple LP. OR course (ineq.form)';
    % Problem formulated on standard form
    A   = [1 2
        4 1 ];
    b_U = [6  12]';
    c   = [-7  -5 ]';
    n   = length(c);
    b_L =-Inf*ones(n,1);
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 8
    Name='Kuhn example.';
    A   = [-2 -9 1 9 1 0
        1/3 1 -1/3 -2 0 1];
    b_U = [0 0 ]';
    c   = [-2 -3 1 12  0 0]';
    %B   = ~[0 0 0 0 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 9
    % Bazaraa page 165
    Name='Beales cycling example.';
    A = [eye(3), [0.25 -8 -1 9; 0.5 -12 -0.5 3; 0 0 1 0]];
    b_U = [0 0 1]';
    c   = [0 0 0 -0.75 20 -0.5 6]';
    %B   = ~[1 1 1 0 0 0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 10
    Name='Winston Ch.6 Review 17. Dual feasible';
    % Could be solved by dual LP solver
    A   = [1  1 -1  0
        1 -2  0 -1 ];
    b_U = [5 8]';
    c   = [2 1 0 0]';
    %B= ~[1 1 0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 11
    %Name='Fletcher 8.21. Redundancy in constraints.';
    Name='Fletcher 8.21. Redundancy in constr.';
    A   = [1  0  2 1 0
        0  1 -1 0 1
        1  1  1 0 0 ];
    b_U = [1  1 2 ]';
    c   = [-1 -1 0 0 0]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 12
    Name='Winston Ex. 4.12 B4. Max | |. Rewritten';
    A   = [4    1    0 0 1 0
        2   -1    0 0 0 1
        2   -3   -1 1 0 0
        ];
    b_U = [4  0.5 0]';
    c   = [0 0  -1 -1 0 0]';
    %B   =~ [0 0 0 1 1 1]';
    n   = length(c);
    b_L = b_U;
    x_L = zeros(n,1);
    x_U = Inf*ones(n,1);
    f_Low=-100;
    x_min=zeros(n,1);
    x_max=5*ones(n,1);
elseif P == 13
    Name='adlittle';
    load lp_probmat P1;
    A   = P1.A;
    b_L = P1.b_L;
    b_U = P1.b_U;
    c   = P1.c;
    x_L = P1.x_L;
    x_U = P1.x_U;
elseif P == 14
    Name='afiro';
    load lp_probmat P2;
    A   = P2.A;
    b_L = P2.b_L;
    b_U = P2.b_U;
    c   = P2.c;
    x_L = P2.x_L;
    x_U = P2.x_U;
elseif P == 15
    Name='bandm';
    load lp_probmat P3;
    A   = P3.A;
    b_L = P3.b_L;
    b_U = P3.b_U;
    c   = P3.c;
    x_L = P3.x_L;
    x_U = P3.x_U;
elseif P == 16
    Name='boeing2';
    load lp_probmat P4;
    A   = P4.A;
    b_L = P4.b_L;
    b_U = P4.b_U;
    c   = P4.c;
    x_L = P4.x_L;
    x_U = P4.x_U;
elseif P == 17
    Name='bore3d';
    load lp_probmat P5;
    A   = P5.A;
    b_L = P5.b_L;
    b_U = P5.b_U;
    c   = P5.c;
    x_L = P5.x_L;
    x_U = P5.x_U;
elseif P == 18
    Name='brandy';
    load lp_probmat P6;
    A   = P6.A;
    b_L = P6.b_L;
    b_U = P6.b_U;
    c   = P6.c;
    x_L = P6.x_L;
    x_U = P6.x_U;
elseif P == 19
    Name='capri';
    load lp_probmat P7;
    A   = P7.A;
    b_L = P7.b_L;
    b_U = P7.b_U;
    c   = P7.c;
    x_L = P7.x_L;
    x_U = P7.x_U;
elseif P == 20
    Name='e226';
    load lp_probmat P8;
    A   = P8.A;
    b_L = P8.b_L;
    b_U = P8.b_U;
    c   = P8.c;
    x_L = P8.x_L;
    x_U = P8.x_U;
elseif P == 21
    Name='etamacro';
    load lp_probmat P9;
    A   = P9.A;
    b_L = P9.b_L;
    b_U = P9.b_U;
    c   = P9.c;
    x_L = P9.x_L;
    x_U = P9.x_U;
elseif P == 22
    Name='farm';
    load lp_probmat P10;
    A   = P10.A;
    b_L = P10.b_L;
    b_U = P10.b_U;
    c   = P10.c;
    x_L = P10.x_L;
    x_U = P10.x_U;
elseif P == 23
    Name='gams10a';
    load lp_probmat P11;
    A   = P11.A;
    b_L = P11.b_L;
    b_U = P11.b_U;
    c   = P11.c;
    x_L = P11.x_L;
    x_U = P11.x_U;
elseif P == 24
    Name='gams30a';
    load lp_probmat P12;
    A   = P12.A;
    b_L = P12.b_L;
    b_U = P12.b_U;
    c   = P12.c;
    x_L = P12.x_L;
    x_U = P12.x_U;
elseif P == 25
    Name='grow7';
    load lp_probmat P13;
    A   = P13.A;
    b_L = P13.b_L;
    b_U = P13.b_U;
    c   = P13.c;
    x_L = P13.x_L;
    x_U = P13.x_U;
elseif P == 26
    Name='israel';
    load lp_probmat P14;
    A   = P14.A;
    b_L = P14.b_L;
    b_U = P14.b_U;
    c   = P14.c;
    x_L = P14.x_L;
    x_U = P14.x_U;
elseif P == 27
    Name='kb2';
    load lp_probmat P15;
    A   = P15.A;
    b_L = P15.b_L;
    b_U = P15.b_U;
    c   = P15.c;
    x_L = P15.x_L;
    x_U = P15.x_U;
elseif P == 28
    Name='kleemin8';
    load lp_probmat P16;
    A   = P16.A;
    b_L = P16.b_L;
    b_U = P16.b_U;
    c   = P16.c;
    x_L = P16.x_L;
    x_U = P16.x_U;
elseif P == 29
    Name='lotfi';
    load lp_probmat P17;
    A   = P17.A;
    b_L = P17.b_L;
    b_U = P17.b_U;
    c   = P17.c;
    x_L = P17.x_L;
    x_U = P17.x_U;

elseif P == 30
    Name='multi';
    load lp_probmat P18;
    A   = P18.A;
    b_L = P18.b_L;
    b_U = P18.b_U;
    c   = P18.c;
    x_L = P18.x_L;
    x_U = P18.x_U;
elseif P == 31
    Name='nug05';
    load lp_probmat P19;
    A   = P19.A;
    b_L = P19.b_L;
    b_U = P19.b_U;
    c   = P19.c;
    x_L = P19.x_L;
    x_U = P19.x_U;
elseif P == 32
    Name='orswq2';
    load lp_probmat P20;
    A   = P20.A;
    b_L = P20.b_L;
    b_U = P20.b_U;
    c   = P20.c;
    x_L = P20.x_L;
    x_U = P20.x_U;
elseif P == 33
    Name='p0033';
    load lp_probmat P21;
    A   = P21.A;
    b_L = P21.b_L;
    b_U = P21.b_U;
    c   = P21.c;
    x_L = P21.x_L;
    x_U = P21.x_U;
elseif P == 34
    Name='p0040';
    load lp_probmat P22;
    A   = P22.A;
    b_L = P22.b_L;
    b_U = P22.b_U;
    c   = P22.c;
    x_L = P22.x_L;
    x_U = P22.x_U;
elseif P == 35
    Name='p0201';
    load lp_probmat P23;
    A   = P23.A;
    b_L = P23.b_L;
    b_U = P23.b_U;
    c   = P23.c;
    x_L = P23.x_L;
    x_U = P23.x_U;
elseif P == 36
    Name='p0282';
    load lp_probmat P24;
    A   = P24.A;
    b_L = P24.b_L;
    b_U = P24.b_U;
    c   = P24.c;
    x_L = P24.x_L;
    x_U = P24.x_U;
elseif P == 37
    Name='p0291';
    load lp_probmat P25;
    A   = P25.A;
    b_L = P25.b_L;
    b_U = P25.b_U;
    c   = P25.c;
    x_L = P25.x_L;
    x_U = P25.x_U;
elseif P == 38
    Name='p0548';
    load lp_probmat P26;
    A   = P26.A;
    b_L = P26.b_L;
    b_U = P26.b_U;
    c   = P26.c;
    x_L = P26.x_L;
    x_U = P26.x_U;
elseif P == 39
    Name='recipe';
    load lp_probmat P27;
    A   = P27.A;
    b_L = P27.b_L;
    b_U = P27.b_U;
    c   = P27.c;
    x_L = P27.x_L;
    x_U = P27.x_U;
elseif P == 40
    Name='refine';
    load lp_probmat P28;
    A   = P28.A;
    b_L = P28.b_L;
    b_U = P28.b_U;
    c   = P28.c;
    x_L = P28.x_L;
    x_U = P28.x_U;
elseif P == 41
    Name='sc105';
    load lp_probmat P29;
    A   = P29.A;
    b_L = P29.b_L;
    b_U = P29.b_U;
    c   = P29.c;
    x_L = P29.x_L;
    x_U = P29.x_U;
elseif P == 42
    Name='sc205';
    load lp_probmat P30;
    A   = P30.A;
    b_L = P30.b_L;
    b_U = P30.b_U;
    c   = P30.c;
    x_L = P30.x_L;
    x_U = P30.x_U;
elseif P == 43
    Name='sc50a';
    load lp_probmat P31;
    A   = P31.A;
    b_L = P31.b_L;
    b_U = P31.b_U;
    c   = P31.c;
    x_L = P31.x_L;
    x_U = P31.x_U;
elseif P == 44
    Name='sc50b';
    load lp_probmat P32;
    A   = P32.A;
    b_L = P32.b_L;
    b_U = P32.b_U;
    c   = P32.c;
    x_L = P32.x_L;
    x_U = P32.x_U;
elseif P == 45
    Name='scagr25';
    load lp_probmat P33;
    A   = P33.A;
    b_L = P33.b_L;
    b_U = P33.b_U;
    c   = P33.c;
    x_L = P33.x_L;
    x_U = P33.x_U;
elseif P == 46
    Name='scagr7';
    load lp_probmat P34;
    A   = P34.A;
    b_L = P34.b_L;
    b_U = P34.b_U;
    c   = P34.c;
    x_L = P34.x_L;
    x_U = P34.x_U;
elseif P == 47
    Name='scfxm1';
    load lp_probmat P35;
    A   = P35.A;
    b_L = P35.b_L;
    b_U = P35.b_U;
    c   = P35.c;
    x_L = P35.x_L;
    x_U = P35.x_U;
elseif P == 48
    Name='scorpion';
    load lp_probmat P36;
    A   = P36.A;
    b_L = P36.b_L;
    b_U = P36.b_U;
    c   = P36.c;
    x_L = P36.x_L;
    x_U = P36.x_U;
elseif P == 49
    Name='sctap1';
    load lp_probmat P37;
    A   = P37.A;
    b_L = P37.b_L;
    b_U = P37.b_U;
    c   = P37.c;
    x_L = P37.x_L;
    x_U = P37.x_U;
elseif P == 50
    Name='share1b';
    load lp_probmat P38;
    A   = P38.A;
    b_L = P38.b_L;
    b_U = P38.b_U;
    c   = P38.c;
    x_L = P38.x_L;
    x_U = P38.x_U;
elseif P == 51
    Name='share2b';
    load lp_probmat P39;
    A   = P39.A;
    b_L = P39.b_L;
    b_U = P39.b_U;
    c   = P39.c;
    x_L = P39.x_L;
    x_U = P39.x_U;
elseif P == 52
    Name='stocfor1';
    load lp_probmat P40;
    A   = P40.A;
    b_L = P40.b_L;
    b_U = P40.b_U;
    c   = P40.c;
    x_L = P40.x_L;
    x_U = P40.x_U;
elseif P == 53
    Name='vtp-base';
    load lp_probmat P41;
    A   = P41.A;
    b_L = P41.b_L;
    b_U = P41.b_U;
    c   = P41.c;
    x_L = P41.x_L;
    x_U = P41.x_U;
elseif P == 54
    Name='zed';
    load lp_probmat P42;
    A   = P42.A;
    b_L = P42.b_L;
    b_U = P42.b_U;
    c   = P42.c;
    x_L = P42.x_L;
    x_U = P42.x_U;
else
    error('lp_prob: Illegal problem number');
end

Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
    [], [], f_Low, x_min, x_max, f_opt, x_opt);
Prob.P = P;

% MODIFICATION LOG:
%
% 981111  hkh  Add computation of x_0 and fix all starting values
%              Corrected error in Beales cycling example
% 981116  hkh  Test on x_0 already computed, to avoid computing twice
% 981208  hkh  Avoid output about x_0, when GUI is initializing
% 990623  hkh  Use qpVarDef and qpProbSet instead of ProbVarDef and ProbSet
% 011030  hkh  Change problem names, avoid too long ones
% 041114  med  Added 42 LP problems
% 041117  med  xxx_prob removed and code added
% 080603  med  Switched to lpAssign, cleaned