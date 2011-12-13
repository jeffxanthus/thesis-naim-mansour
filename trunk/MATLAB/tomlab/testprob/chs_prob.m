% chs_prob: Defines constrained nonlinear problems from Hoch-Schittkowski set
%
% function [probList, Prob] = chs_prob(P);
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

function [probList, Prob] = chs_prob(P,varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'HS 6'...
    ,'HS 7'...
    ,'HS 8'...
    ,'HS 9'...
    ,'HS 10'...
    ,'HS 11'...
    ,'HS 12'...
    ,'HS 13'...
    ,'HS 14'...
    ,'HS 15'...
    ,'HS 16'...
    ,'HS 17'...
    ,'HS 18'...
    ,'HS 19'...
    ,'HS 20'...
    ,'HS 22'...
    ,'HS 23'...
    ,'HS 24'...
    ,'HS 26'...
    ,'HS 27'...
    ,'HS 29'...
    ,'HS 30'...
    ,'HS 31'...
    ,'HS 32'...
    ,'HS 33'...
    ,'HS 34'...
    ,'HS 36'...
    ,'HS 37'...
    ,'HS 39'...
    ,'HS 40'...
    ,'HS 41'...
    ,'HS 42'...
    ,'HS 43'...
    ,'HS 44'...
    ,'HS 46'...
    ,'HS 47'...
    ,'HS 54'...
    ,'HS 55'...
    ,'HS 56'...
    ,'HS 57'...
    ,'HS 59'...
    ,'HS 60'...
    ,'HS 61'...
    ,'HS 62'...
    ,'HS 63'...
    ,'HS 64'...
    ,'HS 65'...
    ,'HS 66. Eckhardt'...
    ,'HS 67'...
    ,'HS 70'...
    ,'HS 71'...
    ,'HS 72'...
    ,'HS 73'...
    ,'HS 74'...
    ,'HS 75'...
    ,'HS 77'...
    ,'HS 78'...
    ,'HS 79'...
    ,'HS 80'...
    ,'HS 81'...
    ,'HS 83'...
    ,'HS 84'...
    ,'HS 86'...
    ,'HS 87'...
    ,'HS 93'...
    ,'HS 95'...
    ,'HS 96'...
    ,'HS 97'...
    ,'HS 98'...
    ,'HS 99'...
    ,'HS 100'...
    ,'HS 101'...
    ,'HS 102'...
    ,'HS 103'...
    ,'HS 104'...
    ,'HS 106'...
    ,'HS 107'...
    ,'HS 108'...
    ,'HS 109'...
    ,'HS 111'...
    ,'HS 112'...
    ,'HS 113'...
    ,'HS 114'...
    ,'HS 116'...
    ,'HS 117'...
    ,'HS 119'...
    ,'HS 215'...
    ,'HS 216'...
    ,'HS 217'...
    ,'HS 218'...
    ,'HS 219'...
    ,'HS 220'...
    ,'HS 221'...
    ,'HS 222'...
    ,'HS 223'...
    ,'HS 225'...
    ,'HS 226'...
    ,'HS 227'...
    ,'HS 228'...
    ,'HS 230'...
    ,'HS 231'...
    ,'HS 232'...
    ,'HS 233'...
    ,'HS 234'...
    ,'HS 235'...
    ,'HS 236'...
    ,'HS 237'...
    ,'HS 238'...
    ,'HS 239'...
    ,'HS 248'...
    ,'HS 249'...
    ,'HS 250'...
    ,'HS 251'...
    ,'HS 252'...
    ,'HS 253'...
    ,'HS 254'...
    ,'HS 263'...
    ,'HS 264'...
    ,'HS 265'...
    ,'HS 270'...
    ,'HS 284'...
    ,'HS 285'...
    ,'HS 315'...
    ,'HS 316'...
    ,'HS 317'...
    ,'HS 318'...
    ,'HS 319'...
    ,'HS 320'...
    ,'HS 321'...
    ,'HS 322'...
    ,'HS 323'...
    ,'HS 324'...
    ,'HS 325'...
    ,'HS 326'...
    ,'HS 327'...
    ,'HS 329'...
    ,'HS 330'...
    ,'HS 331'...
    ,'HS 332'...
    ,'HS 335'...
    ,'HS 336'...
    ,'HS 337'...
    ,'HS 338'...
    ,'HS 339'...
    ,'HS 340'...
    ,'HS 341'...
    ,'HS 342'...
    ,'HS 343'...
    ,'HS 344'...
    ,'HS 345'...
    ,'HS 346'...
    ,'HS 347'...
    ,'HS 353'...
    ,'HS 354'...
    ,'HS 355'...
    ,'HS 356'...
    ,'HS 360'...
    ,'HS 361'...
    ,'HS 365'...
    ,'HS 366'...
    ,'HS 367'...
    ,'HS 369'...
    ,'HS 372'...
    ,'HS 373'...
    ,'HS 374'...
    ,'HS 375'...
    ,'HS 376'...
    ,'HS 377'...
    ,'HS 378'...
    ,'HS 380'...
    ,'HS 382'...
    ,'HS 383'...
    ,'HS 384'...
    ,'HS 385'...
    ,'HS 386'...
    ,'HS 387'...
    ,'HS 388'...
    ,'HS 389'...
    ,'HS 394'...
    ,'HS 395'...
    ); % CHANGE: MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

if P < 46
    Prob = chs_prob1(P);
elseif P < 91
    Prob = chs_prob2(P);
elseif P < 136
    Prob = chs_prob3(P);
elseif P < 181
    Prob = chs_prob4(P);
else
    error('chs_prob: Illegal problem number');
end

function Prob = chs_prob1(P)

HessPattern = []; pSepFunc = []; ConsPattern = []; f_Low = [];

if P == 1
    Name='HS 6';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [-1.2;1];
    x_L = []; 
    x_U = [];
    f_Low = 0;
    x_min = [-2 -2];
    x_max = [ 2  2];
    x_opt = [1 1];
    f_opt = 0;
elseif P==2
    Name='HS 7';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2];
    x_max = [ 2  2];
    x_opt = [0 sqrt(3)];
    f_opt = -sqrt(3);
elseif P==3
    Name='HS 8';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [2;1];
    x_L = []; 
    x_U = [];
%     x_min = [-2 -2];
%     x_max = [ 2  2];
    x_min = [-5 -5];
    x_max = [ 5  5];
    %x_opt = [];
    % 4 optimal solutions, intersections of 2 constraints
    x1 = sqrt(25/2 + sqrt(301)/2);
    x2 = 9/x1;
    x_opt = [x1 x2 ; x2 x1 ; -x1 -x2 ; -x2 -x1];
    f_opt = -1;
elseif P==4
    Name='HS 9';
    b_L = 0;
    b_U = 0;
    A   = [4 -3]; 
    c_L = []; c_U = [];
    x_0 = [0;0];
    x_L = []; 
    x_U = [];
%     x_min = [-2 -2];
%     x_max = [ 2  2];
    x_min = [-10000 -10000];
    x_max = [0 0];
    %x_opt = [];
    % x_opt found using glcCluster. The objective function
    % involves sin and cos, probably multiple optimas.
    x_opt = [-6315 -8420 ; -5619 -7492];
    f_opt = -0.5;
elseif P==5
    Name='HS 10';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [-10;10];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2];
    x_max = [ 2  2];
    x_opt = [0 1];
    f_opt = -1;
elseif P==6
    Name='HS 11';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [4.9;.1];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2];
    x_max = [ 2  2];
    %a = 7.5*sqrt(6)+sqrt(338.5);
    % This is not correct. f(x_opt) >>> f_opt.
    % f_opt seems to be correct, so probably something
    % wrong with the analytical expression for x_opt.
    %x_opt = [(a-1/a)/sqrt(6) (a^2-2+1/a^2)/6];
    %f_opt = -8.498464223;
    
    % Solved with glcCluster:
    x_opt = [1.234772892632396 1.524664093087327];
    f_opt = -8.498464233194383;
elseif P==7
    Name='HS 12';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [0;0];
    x_L = []; 
    x_U = [];
%     x_min = [-2 -2];
%     x_max = [ 2  2];
    x_min = [-2 -2];
    x_max = [ 5  5];
    x_opt = [2 3];
    f_opt = -30;
elseif P==8
    Name='HS 13';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [-2;-2];
    x_L = [0;0]; 
    x_U = [];
    x_min = [-4 -4];
    x_max = [ 4  4];
    x_opt = [1 0];
    f_opt = 1;
elseif P==9
    % Schittkowski 14. Bracken, McCormick, Himmelblau. Start (2,2), f=1
    Name='HS 14';
    A=[1 -2];
    b_L=-1;
    b_U=-1;
    c_L=0;
    c_U=[];
    
    x_0=[2;2];  % Not feasible
    x_opt=[0.5*(sqrt(7)-1) 0.25*(sqrt(7)+1)]; % 0.822875656 0.911437828
    f_opt=9-2.875*sqrt(7); % 1.3934649806893
    x_L=[-5;-5];
    x_U=[5;5];
    f_Low=0;
    x_min=[0;0];
    x_max=[2;2];
elseif P==10
    Name='HS 15';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [-2;1];
    x_L = []; 
    x_U = [.5;inf];
    x_min = [-4 -4];
    x_max = [ 4  4];
    x_opt = [.5 2];
    f_opt = 306.5;
elseif P==11
    Name='HS 16';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [-2;1];
    x_L = [-.5;-inf]; 
    x_U = [.5;1];
    x_min = [-3 -3];
    x_max = [ 3  3];
    x_opt = [.5 .25];
    f_opt = .25;
elseif P==12
    Name='HS 17';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [-2;1];
    x_L = [-.5;-inf]; 
    x_U = [.5;1];
    x_min = [-4 -4];
    x_max = [ 4  4];
    x_opt = [0 0];
    f_opt = 1;
elseif P==13
    Name='HS 18';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [2;2];
    x_L = [2;0]; 
    x_U = [50;50];
    x_min = [-4 -4];
    x_max = [ 4  4];
    x_opt = [sqrt(250) sqrt(2.5)];
    f_opt = 5;
elseif P==14
    Name='HS 19';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [20.1;5.84];
    x_L = [13;0]; 
    x_U = [100;100];
    x_min = [0 0];
    x_max = [ 25  10];
    x_opt = [14.095 .84296079];
    f_opt = -6961.81381;
elseif P==15
    Name='HS 20';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0;0];
    c_U = [];
    x_0 = [-2;1];
    x_L = [-.5;-inf]; 
    x_U = [.5;inf];
    x_min = [-2 -2];
    x_max = [ 2  2];
    x_opt = [.5 .5*sqrt(3)];
    f_opt = 81.5-25*sqrt(3);
elseif P==16
    Name='HS 22';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [2;2];
    x_L = []; 
    x_U = [];
    x_min = [0 0];
    x_max = [ 3  3];
    x_opt = [1 1];
    f_opt = 1;
elseif P==17
    Name='HS 23';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0;0;0;0];
    c_U = [];
    x_0 = [3;1];
    x_L = [-50;-50]; 
    x_U = [50;50];
    x_min = [0 0];
    x_max = [ 4  4];
    x_opt = [1 1];
    f_opt = 2;
elseif P==18
    % Schittkowski 24. Betts and Box. Start in (1,.5), f= -.01336459
    Name='HS 24';
    % 3 linear inequalities
    r=sqrt(3);
    A=[1/r -1; 1 r; -1 -r];
    b_L=[0 0 -6]';
    b_U=inf*ones(3,1);
    c_L=[]; c_U=[];
    x_0=[1;0.5]; % Feasible
    x_L=[0;0];
    x_U=[5;5];
    f_Low=-10;
    x_min=[0;0];
    x_max=[4;2];
    x_opt=[3 sqrt(3)];
    f_opt=-1;
elseif P==19
    Name='HS 26';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [-2.6;2;2];
    x_L = []; 
    x_U = [];
    x_min = [-3 -3 -3];
    x_max = [ 3 3 3];
    % Something is wrong here!!
    % If c_L = c_U = -2, then x_opt = [1 1 1] is correct.
    %x_opt = [1 1 1];
    
    % If c_L = c_U = 0, then this is correct:
    xx1 =  1.475096131425521;
    xx2 = -1.237247540050211;
    x_opt = [xx1 ; xx2]*ones(1,3);
    f_opt = 0;
elseif P==20
    Name='HS 27';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2 2 2];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2];
    x_max = [ 2 2 2];
    x_opt = [-1 1 0];
    f_opt = 0.04;
elseif P==21
    Name='HS 29';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1;1;1];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0];
    x_max = [4 4 4];
    x_opt = [4 2*sqrt(2) 2];
    f_opt = -16*sqrt(2);
elseif P==22
    Name='HS 30';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1;1;1];
    x_L = [1;-10;-10]; 
    x_U = [10;10;10];
    x_min = [0 0 0];
    x_max = [2 2 2];
    x_opt = [1 0 0];
    f_opt = 1;
elseif P==23
    Name='HS 31';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1;1;1];
    x_L = [-10;1;-10]; 
    x_U = [10;10;1];
    x_min = [-1 -1 -1];
    x_max = [1 1 1];
    x_opt = [1/sqrt(3) sqrt(3) 0];
    f_opt = 6;
elseif P==24
    Name='HS 32';
    b_L=1; b_U=1; A=[1 1 1]; 
    c_L = 0;
    c_U = [];
    x_0 = [.1;.7;.2];
    x_L = [0;0;0]; 
    x_U = [];
    x_min = [-1 -1 -1];
    x_max = [2;2;2];
    x_opt = [0 0 1];
    f_opt = 1;
elseif P==25
    Name='HS 33';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0;0;3];
    x_L = [0;0;0]; 
    x_U = [inf;inf;5];
    x_min = [-1 -1 -1];
    x_max = [2;2;2];
    x_opt = [0 sqrt(2) sqrt(2)];
    f_opt = sqrt(2)-6;
elseif P==26
    Name='HS 34';
    b_L=[]; b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0;1.05;2.9];
    x_L = [0;0;0]; 
    x_U = [100;100;10];
    x_min = [-1 -1 -1];
    x_max = [2;2;2];
    x_opt = [log(log(10)) log(10) 10];
    f_opt = -log(log(10));
elseif P==27
    Name='HS 36';
    b_L=[]; b_U=[]; A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [10;10;10];
    x_L = [0;0;0]; 
    x_U = [20;11;42];
    x_min = [5 5 5];
    x_max = [20;20;20];
    x_opt = [20 11 15];
    f_opt = -3300;
elseif P==28
    Name='HS 37';
    b_L=0; b_U=72; A=[1 2 2]; 
    c_L = [];
    c_U = [];
    x_0 = [10;10;10];
    x_L = [0;0;0]; 
    x_U = [42;42;42];
    x_min = [5 5 5];
    x_max = [25;25;25];
    x_opt = [24 12 12];
    f_opt = -3456;
elseif P==29
    Name='HS 39';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [2;2;2;2];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0];
    x_max = [2 2 2 2];
    x_opt = [1 1 0 0];
    f_opt = -1;
elseif P==30
    Name='HS 40';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [.8;.8;.8;.8];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0];
    x_max = [2 2 2 2];
    % Sign error!! x_opt(3) must also be positive
    x_opt = [2^(-1/3) 2^(-1/2) 2^(-11/12) 2^(-1/4)];
    f_opt = -.25;
elseif P==31
    Name='HS 41';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2;2;2];
    x_L = [0;0;0;0]; 
    x_U = [1 1 1 2];
    x_min = [0 0 0 0];
    x_max = [2 2 2 2];
    x_opt = [2/3 1/3 1/3 2];
    f_opt = 52/27;
elseif P==32
    Name='HS 42';
    b_L=2;b_U=2;A=[1 0 0 0]; 
    c_L = 0;
    c_U = 0;
    x_0 = [1;1;1;1];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0];
    x_max = [2 2 2 2];
    % Small mistake, x_opt(3) should be divided by 10
    %x_opt = [2 2 6*sqrt(2) .8*sqrt(2)];
    x_opt = [2 2 .6*sqrt(2) .8*sqrt(2)];
    f_opt = 28-10*sqrt(2);
elseif P==33
    Name='HS 43';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [];
    x_0 = [0;0;0;0];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2 -2];
    x_max = [ 2  2  2  2];
    x_opt = [0 1 2 -1];
    f_opt = -44;
elseif P==34
    Name='HS 44';
    b_L=-inf*ones(6,1);
    b_U=[8;12;12;8;8;5];
    A=[1 2 0 0
        4 1 0 0
        3 4 0 0
        0 0 2 1
        0 0 1 2
        0 0 1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [0;0;0;0];
    x_L = [0;0;0;0]; 
    x_U = [];
    x_min = [-1 -1 -1 -1];
    x_max = [4 4 4 4];
    x_opt = [0 3 0 4];
    f_opt = -15;
elseif P==35
    Name='HS 46';
    b_L=[];
    b_U=[];
    A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [.5*sqrt(2);1.75;.5;2;2];
    x_L = []; 
    x_U = [];
    x_min = [-1 -1 -1 -1 -1];
    x_max = [4 4 4 4 4];
    x_opt = [1 1 1 1 1];
    f_opt = 0;
elseif P==36
    Name='HS 47';
    b_L=[];
    b_U=[];
    A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [2;sqrt(2);-1;2-sqrt(2);.5];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0 0];
    x_max = [4 4 4 4 4];
    % Improved optima using glcCluster
%     x_opt = [1 1 1 1 1];
%     f_opt = 0;
    x_opt = [0.677004616701280 0.726089640962586 ...
             1.215491109042481 1.751328995197990 1.477094801614824];
    f_opt = -0.026714182694062;
elseif P==37
    Name='HS 54';
    b_L=1.76e4;
    b_U=1.76e4;
    A=[1 4e3 0 0 0 0]; 
    c_L = [];
    c_U = [];
    x_0 = [6e3;1.5;4e6;2;3e-3;5e7];
    x_L = [0;-10;0;0;-1;0]; 
    x_U = [2e4;10;1e7;20;1;2e8];
    x_min = [0 0 0 0 0 0];
    x_max = [4 4 4 4 4 4];
    % Something is wrong!!  chs_f(x_opt) ~= f_opt
    % x_opt fulfills linear constraints. Maybe f_opt is wrong??
    % chs_f(x_opt) --> f_opt = 0.
    x_opt = [91600/7 79/70 2e6 10 1e-3 1e8];
    f_opt = 0;
    %f_opt = -exp(-27/280);
elseif P==38
    Name='HS 55';
    b_L=[6 3 2 1 2 2];
    b_U=[6 3 2 1 2 2];
    A=[1 2 0 0 5 0
        1 1 1 0 0 0
        0 0 0 1 1 1
        1 0 0 1 0 0
        0 1 0 0 1 0
        0 0 1 0 0 1]; 
    c_L = [];
    c_U = [];
    x_0 = [1;2;0;0;0;2];
    x_L = [0;0;0;0;0;0]; 
    x_U = [1;inf;inf;1;inf;inf];
    x_min = [0 0 0 0 0 0];
    x_max = [4 4 4 4 4 4];
    x_opt = [0 4/3 5/3 1 2/3 1/3];
    f_opt = 19/3;
elseif P==39
    Name='HS 56';
    b_L=[];
    b_U=[];
    A=[]; 
    c_L = [0;0;0;0];
    c_U = [0;0;0;0];
    a=asin(sqrt(1/4.2));
    b=asin(sqrt(5/7.2));
    c=asin(sqrt(4/7));
    d=asin(sqrt(2/7));
    x_0 = [1;1;1;a;a;a;b];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0 0 0 0];
    x_max = [4 4 4 4 4 4 4];
    % Possible multiple minima
    x_opt = [2.4 1.2 1.2 c d d .5*pi];
    f_opt = -3.456;
elseif P==40
    Name='HS 57';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [.42;5];
%     x_L = [.4;-4];
%     x_U = [];
%     x_min = [0 0];
%     x_max = [1 1];
    x_L = [];
    x_U = [];
    x_min = [0 0];
    x_max = [2 2];
    x_opt = [.419952675 1.284845629];
    f_opt = .02845966972;
elseif P==41
    Name='HS 59';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [90;10];
    x_L = [0;0]; 
    x_U = [75;65];
    x_min = [0 0];
    x_max = [10 10];
%     x_opt = [13.55010424 51.66018129];
%     f_opt = -7.804226324;
    x_opt = [13.551245575708313  51.655785872749753];
    f_opt = -7.802785527345220;
elseif P==42
    Name='HS 60';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2;2];
    x_L = [-10;-10;-10]; 
    x_U = [10;10;10];
    x_min = [0 0 0];
    x_max = [4 4 4];
    x_opt = [1.104859024 1.196674194 1.535262257];
    f_opt = .03256820025;
elseif P==43
    Name='HS 61';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [0;0;0];
    x_L = []; 
    x_U = [];
    x_min = [0 -3 0];
    x_max = [6  6 6];
    x_opt = [5.326770157 -2.118998639 3.210464239];
    f_opt = -143.6461422;
elseif P==44
    Name='HS 62';
    b_L=1;b_U=1;A=[1 1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [.7;.2;.1];
    x_L = [0;0;0]; 
    x_U = [1;1;1];
    x_min = [0 0 0];
    x_max = [1 1 1];
    x_opt = [.6178126908 .328202223 .05398508606];
    f_opt = -26272.51448;
elseif P==45
    Name='HS 63';
    b_L=56;b_U=56;A=[8 14 7]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2;2];
    x_L = [0;0;0]; 
    x_U = [inf;inf;inf];
    x_min = [0 0 0];
    x_max = [4 4 4];
    x_opt = [3.512118414 .2169881741 3.552174034];
    f_opt = 961.7151721;
end

c  = 'chs_c';
dc = 'chs_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('chs_f','chs_g','chs_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P      = P;

function Prob = chs_prob2(P)

HessPattern = []; pSepFunc = []; ConsPattern = []; f_Low = [];
f_opt = [];

if P==46
    Name='HS 64';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1;1;1];
    x_L = [1e-5;1e-5;1e-5]; 
    x_U = [inf;inf;inf];
    x_min = [0 0 0];
    x_max = [110 100 210];
    x_opt = [108.7347175;85.12613942;204.3247078]';
    f_opt = 6299.842428;
elseif P==47
    Name='HS 65';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [-5;5;0];
    x_L = [-4.5;-4.5;-5]; 
    x_U = [4.5;4.5;5];
    x_min = [-5 -5 -5];
    x_max = [5 5 5];
    x_opt = [3.650461821;3.65046168;4.6204170507]';
    f_opt = .9535288567;
elseif P==48
    % HS 66. Eckhardt. Start in (0,1.05,2.9). f=.58
    Name='HS 66. Eckhardt';
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
    x_opt=[0.1841264879; 1.202167873; 3.327322322]';
    f_opt=.5181632741;
elseif P==49
    Name='HS 67';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [1745;12000;110];
    x_L = [1e-5;1e-5;1e-5]; 
    x_U = [2e3;1.6e4;1.2e2];
    x_min = [-5 -5 -5];
    x_max = [5 5 5];
    x_opt = [1728.371286;16000.000000;98.14151402]';
    f_opt = -1162.036507;
elseif P==50
    Name='HS 70';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [2;4;.04;2];
    x_L = [.00001;.00001;.00001;.00001]; 
    x_U = [100;100;1;100];
    x_min = [-5 -5 -5 -5];
    x_max = [5 5 5 5];
    x_opt = [12.27695;4.631788;.3128625;2.029290]';
    f_opt = .007498464;
elseif P==51
    Name='HS 71';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;0];
    x_0 = [1;5;5;1];
    x_L = [1;1;1;1]; 
    x_U = [5;5;5;5];
    x_min = [-5 -5 -5 -5];
    x_max = [5 5 5 5];
    x_opt = [1;4.7429994;3.8211503;1.3794082]';
    f_opt = 17.0140173;
elseif P==52
    Name='HS 72';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [1;1;1;1];
    x_L = [.001;.001;.001;.001]; 
    x_U = [4e5;3e5;2e5;1e5];
    x_min = [-5 -5 -5 -5];
    x_max = [5 5 5 5];
    x_opt = [193.4071;179.5475;185.0186;168.7062]';
    f_opt = 727.67937;
elseif P==53
    Name='HS 73';
    b_L=1;b_U=1;A=[1 1 1 1]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [1;1;1;1];
    x_L = [0;0;0;0]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-1 -1 -1 -1];
    x_max = [1 1 1 1];
    x_opt = [.6355216;-.12e-11;.3127019;.05177655]';
    f_opt = 29.894378;
elseif P==54
    Name='HS 74';
    b_L=[-.55;-.55];b_U=[inf;inf];A=[0 0 -1 1;0 0 1 -1]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [0;0;0;0];
    x_L = [0;0;-.55;-.55]; 
    x_U = [1200;1200;.55;.55];
    x_min = [-1 -1 -1 -1];
    x_max = [1 1 1 1];
    x_opt = [679.9453;1026.067;.110879;-.3962336]';
    f_opt = 5126.4981;
elseif P==55
    Name='HS 75';
    b_L=[-.48;-.48];b_U=[inf;inf];A=[0 0 -1 1;0 0 1 -1]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [0;0;0;0];
    x_L = [0;0;-.48;-.48]; 
    x_U = [1200;1200;.48;.48];
    x_min = [-1 -1 -1 -1];
    x_max = [1 1 1 1];
    x_opt = [776.1592;925.1949;.05110879;-.4288911]';
    f_opt = 5174.4129;
elseif P==56
    Name='HS 77';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [2;2;2;2;2];
    x_L = []; 
    x_U = [];
    x_min = [-1 -1 -1 -1 -1];
    x_max = [ 2  2  2  2  2];
    x_opt = [1.166172;1.182111;1.380257;1.506036;.6109203]';
    f_opt = .24150513;
elseif P==57
    Name='HS 78';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [-2;1.5;2;-1;-1];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2];
    x_opt = [-1.717142;1.595708;1.827248;-.7636429;-7636435]';
    f_opt = -2.91970041;
elseif P==58
    Name='HS 79';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [2;2;2;2;2];
    x_L = []; 
    x_U = [];
    x_min = [0 0 0 0 0];
    x_max = [2 2 2 2 2];
    x_opt = [1.191127;1.362603;1.472818;1.635017;1.679081]';
    f_opt = .0787768209;
elseif P==59
    Name='HS 80';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [-2;2;2;-1;-1];
    x_L = [-2.3;-2.3;-3.2;-3.2;-3.2]; 
    x_U = [2.3;2.3;3.2;3.2;3.2];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2];
    x_opt = [-1.717143;1.595709;1.827247;-.7636413;-.7636450]';
    f_opt = .0539498478;
elseif P==60
    Name='HS 81';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [0;0;0];
    x_0 = [-2;2;2;-1;-1];
    x_L = [-2.3;-2.3;-3.2;-3.2;-3.2]; 
    x_U = [2.3;2.3;3.2;3.2;3.2];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2];
    % Something is wrong here!!  chs_c(x_opt) not feasible. 
    % Small typo in x_opt(2). Change from 1.159571 to 1.595710
    % Also found alternative optima.
    x_opt = [-1.717142 1.595710 1.827248 -.7636474 -.7636390 ;
             -1.717142 1.595710 1.827248  .7636474  .7636390];
    f_opt = .0539498478;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif P==61
    Name='HS 83';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [92;20;5];
    x_0 = [78;33;27;27;27];
    x_L = [78;33;27;27;27]; 
    x_U = [102;45;45;45;45];
    x_min = [55 32 26 26 25];
    x_max = [80 35 46 37 35];
    x_opt = [78;33;29.99526;45;36.77581]';
    f_opt = -30665.53867;
    
elseif P==62
    Name='HS 84';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [294000;294000;277200];
    x_0 = [2.52;2;37.5;9.25;6.8];
    x_L = [0;1.2;20;9;6.5]; 
    x_U = [1000;2.4;60;9.3;7];
    x_min = [0 1 10 0 0];%Marginal used only for grafic presentation
    x_max = [150 10 100 15 15];%Marginal used only for grafic presentation
    x_opt = [4.53743097;2.4;60;9.3;7];
    f_opt = -5280335.133;
elseif P==63
    % 10 linear inequalities
    Name='HS 86';
    b_L=[-40 -2  -0.25 -4  -4  -1  -40  -60  5  1];
    b_U=inf*ones(10,1);
    A = [-16   2   0     1  0;
        0  -2   0     4  2;
        -3.5  0   2     0  0;
        0  -2   0    -4 -1;
        0  -9  -2     1 -2.8;
        2   0  -4     0  0;
        -1  -1  -1    -1 -1;
        -1  -2  -3    -2 -1;
        1   2   3     4  5;
        1   1   1     1  1];          
    c_L = [];   
    c_U = [];
    x_0 = [0;0;0;0;1]; %Feasible
    x_L = [0;0;0;0;0];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1];
    x_max = [1;1;1;1;1];
    x_opt = [0.3;0.33346761;0.4;0.42831010;0.22396487]';
    f_opt = -32.34867897;
    
elseif P== 64
    Name='HS 87';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0];
    c_U = [0;0;0;0];
    x_0 = [390;1000;419.5;340.5;198.175;0.5];
    x_L = [0;0;340;340;-1000;0]; 
    x_U = [400;1000;420;420;10000;0.5236];
    x_min = [0;0;0;0;-100;0];%Marginal used only for grafic presentation
    x_max = [150;250;400;450;300;10];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [107.8119;196.3186;373.8307;420;213.0713;0.1532920]';
%     f_opt = 8927.5977;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [  0  311.3851364441459  356.8038659421023 ...
             420  88.8228659135595   0.2570430426130];
    f_opt = 8718.783820436085;
elseif P== 65
    Name='HS 93';
    % 2 equalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [5.54;4.4;12.02;11.82;0.702;0.852];
    x_L = [0;0;0;0;0;0]; 
    x_U = [];
    x_min = [0;0;0;0;0;0];%Marginal used only for grafic presentation
    x_max = [10;10;15;15;10;10];%Marginal used only for grafic presentation
    x_opt = [5.332666;4.656744;10.43299;12.08230;0.7526074;0.87865084]';
    f_opt = 135.075961;
elseif P== 66
    Name='HS 95';
    b_L=[];b_U=[];A=[]; 
    c_L = [4.97;-1.88;-29.08;-78.02];
    c_U = [];
    x_0 = [0;0;0;0;0;0];
    x_L = [0;0;0;0;0;0]; 
    x_U = [0.31;0.46;0.068;0.042;0.028;0.0134];
    x_min = [-5;-5;-5;-5;-5;-5];%Marginal used only for grafic presentation
    x_max = [5;5;5;5;5;5];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [0;0;0;0;0;0]';
%     f_opt = 0.015619514;
    % Used glcCluster to find new feasible x_opt and f_opt.
    x_opt = [0 0 0 0 0 0.003323303243063]';
    f_opt = 0.015619525242394;
    
elseif P== 67
    Name='HS 96';
    b_L=[];b_U=[];A=[]; 
    c_L = [4.97;-1.88;-69.08;-118.02];
    c_U = [];
    x_0 = [0;0;0;0;0;0];
    x_L = [0;0;0;0;0;0]; 
    x_U = [0.31;0.46;0.068;0.042;0.028;0.0134];
    x_min = [-5;-5;-5;-5;-5;-5];%Marginal used only for grafic presentation
    x_max = [5;5;5;5;5;5];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [0;0;0;0;0;0]';
%     f_opt = 0.015619514;
    % Used glcCluster to find new feasible x_opt and f_opt.
    x_opt = [0 0 0 0 0 0.003323303243063]';
    f_opt = 0.015619525242394;
    
elseif P== 68
    Name='HS 97';
    b_L=[];b_U=[];A=[]; 
    c_L = [32.97;25.12;-29.08;-78.02];
    c_U = [];
    x_0 = [0;0;0;0;0;0];
    x_L = [0;0;0;0;0;0]; 
    x_U = [0.31;0.46;0.068;0.042;0.028;0.0134];
    x_min = [-5;-5;-5;-5;-5;-5];%Marginal used only for grafic presentation
    x_max = [5;5;5;5;5;5];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [0;0;0;0;0;0.0033233033]';
%     f_opt = 3.1358091;
    % Used glcCluster to find new feasible x_opt and f_opt.
    x_opt = [0.268564912271566 0 0 0 0.028 0.0134];
    f_opt = 3.135809122767733;
    
elseif P== 69
    Name='HS 98';
    b_L=[];b_U=[];A=[]; 
    c_L = [32.97;25.12;-124.08;-173.02];
    c_U = [];
    x_0 = [0;0;0;0;0;0];
    x_L = [0;0;0;0;0;0]; 
    x_U = [0.31;0.46;0.068;0.042;0.028;0.0134];
    x_min = [-5;-5;-5;-5;-5;-5];%Marginal used only for grafic presentation
    x_max = [5;5;5;5;5;5];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [0;0;0;0;0;0.0033233033]';
%     f_opt = 3.1358091;
    % Used glcCluster to find new feasible x_opt and f_opt.
    x_opt = [0.268564912274200 0 0 0 0.028 0.0134];
    f_opt = 3.135809122779061;
elseif P== 70
    Name='HS 99';
    % 2 nonlinear equalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = 0.5*ones(7,1);%% Not feasible
    x_L = zeros(7,1); 
    x_U = 1.58*ones(7,1);
    x_min = -2*ones(7,1);%Marginal used only for grafic presentation
    x_max = 2*ones(7,1);%Marginal used only for grafic presentation
    %x_opt = [0.5424603;0.5290159;0.5084506;0.4802693;0.4512352;0.4091878;0.3527847]';
    %f_opt = -0.831079892*10^9;
    % Used glcCluster to find new, more feasible, x_opt and f_opt.
    x_opt = [0.542467805835978  0.529021417683741  0.508449153351392 ...
             0.480268853512796  0.451236354555062  0.409183086454820 ...
             0.352787887841798];
    f_opt = -8.310798915101076e+008;
elseif P== 71
    Name='HS 100';
    % 4 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0];
    c_U = [];
    x_0 = [1;2;0;4;0;1;1];%% Not feasible
    x_L = []; 
    x_U = [];
    x_min = -4*ones(7,1);
    x_max = 5*ones(7,1);
    x_opt = [2.330499;1.951372;-0.4775414;4.365726;-0.6244870;1.038131;1.594227]';
    f_opt = 680.6300573;
elseif P== 72
    Name='HS 101';
    % 5 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;100];
    c_U = [inf;inf;inf;inf;3000];
    x_0 = [6;6;6;6;6;6;6];%% Not feasible
    x_L = [0.1;0.1;0.1;0.1;0.1;0.1;0.01]; 
    x_U = [10;10;10;10;10;10;10];
    x_min = [-4;-4;-4;-4;-4;-4;-4];%Marginal used only for grafic presentation
    x_max = [4;4;4;4;4;4;4];%Marginal used only for grafic presentation
    x_opt = [2.856159;0.6108230;2.150813;4.712874;0.9994875;1.347508;0.3165277]';
    f_opt = 1809.76476;
elseif P== 73
    Name='HS 102';
    % 5 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;100];
    c_U = [inf;inf;inf;inf;3000];
    x_0 = [6;6;6;6;6;6;6];%% Not feasible
    x_L = [0.1;0.1;0.1;0.1;0.1;0.1;0.01]; 
    x_U = [10;10;10;10;10;10;10];
    x_min = [-4;-4;-4;-4;-4;-4;-4];%Marginal used only for grafic presentation
    x_max = [4;4;4;4;4;4;4];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [2.856159;0.6108230;2.150813;4.712874;0.9994875;1.347508;0.3165277]';
%     f_opt = 991.880571;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [3.896260818291297 0.809358142205370 2.664380149037519 ...
             4.300885927707119 0.853552241316500 1.095283114124308 ...
             0.027310186429933];
    f_opt = 911.88057147062705;
elseif P== 74
    Name='HS 103';
    % 5 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;100];
    c_U = [inf;inf;inf;inf;3000];
    x_0 = [6;6;6;6;6;6;6];%% Not feasible
    x_L = [0.1;0.1;0.1;0.1;0.1;0.1;0.01]; 
    x_U = [10;10;10;10;10;10;10];
    x_min = [-4;-4;-4;-4;-4;-4;-4];%Marginal used only for grafic presentation
    x_max = [4;4;4;4;4;4;4];%Marginal used only for grafic presentation
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [2.856159;0.6108230;2.150813;4.712874;0.9994875;1.347508;0.3165277]';
%     f_opt = 543.667958;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [4.394103578721100 0.854468790267496 2.843229650594898 ...
             3.399979053896162 0.722926036238323 0.870406391771914 ...
             0.024638832659986];
    f_opt = 543.6679569438567;
elseif P== 75
    Name='HS 104';
    % 4 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;1];
    c_U = [inf;inf;inf;inf;4.2];
    x_0 = [6;3;0.4;0.2;6;6;1;0.5];%% Not feasible
    x_L = [0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1]; 
    x_U = [10;10;10;10;10;10;10;10];
    x_min = [-4;-4;-4;-4;-4;-4;-4;-4];%Marginal used only for grafic presentation
    x_max = [10;10;10;10;10;10;10;10];%Marginal used only for grafic presentation
    x_opt = [6.465114;2.232709;0.6673975;0.5957564;5.932676;5.527235;1.013322;0.4006682]';
    f_opt = 3.9511634396;
elseif P== 76
    Name='HS 106';
    % 3 linear inequalities
    % 3 nonlinear inequalities
    b_L=[-1 -1 -1];
    b_U=inf*ones(3,1);
    A=[0 0 0 -0.0025   0    -0.0025     0  0;
        0 0 0  0.0025  -0.0025  0   -0.0025 0;
        0 0 0   0       0.01    0     0  -0.01]; 
    c_L = [0;0;0];
    c_U = [];
    x_0 = [5000;5000;5000;200;350;150;225;425];%% Not feasible
    x_L = [100;1000;1000;10;10;10;10;10]; 
    x_U = [10000;10000;10000;1000;1000;1000;1000;1000];
    x_min = [100;1000;1000;4000;100;100;100;100];
    x_max = [600;2000;4000;6000;400;300;400;500];
    x_opt = [579.3167;1359.943;5110.071;182.0174;295.5985;217.9799;286.4162;395.5979]';
    f_opt = 7049.330923;
elseif P== 77
    Name='HS 107';
    % 6 nonlinear equalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0];
    c_U = [0;0;0;0;0;0];
    x_0 = [0.8;0.8;0.2;0.2;1.0454;1.0454;0;0;0];%Not feasible
    x_L = [0;0;-inf;-inf;0.90909;0.90909;0.90909;-inf;-inf]; 
    x_U = [inf;inf;inf;inf;1.0909;1.0909;1.0909;inf;inf];  
    x_min = [-2;-2;-2;-2;-2;-2;-2;-2;-2];
    x_max = [2;2;2;2;2;2;2;2;2];
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [0.6670095;1.022388;0.2238879;0.1848217;1.090900;1.090900;
%              1.069360;0.1066126;-0.3387867]';
%     f_opt = 5055.011803;
    % Used glcCluster to find new feasible x_opt and f_opt.
    x_opt = [ 0.910936935877514  1.032536048084956 0.784854183570069 ...
              1.730383672818811  0.767150404858614 1.055042057545561 ...
             -0.559437868334162 -0.131798460445879 2];
    f_opt = 6287.663211231529;
elseif P== 78
    Name='HS 108';
    % 13 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [1;1;1;1;1;1;1;1;1];%% Not feasible
    x_L = [-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;0];
    x_U = []; 
    x_min = [-2;-2;-2;-2;-2;-2;-2;-2;-2];
    x_max = [2;2;2;2;2;2;2;2;2];
    x_opt = [0.8841292;0.4672425;0.03742076;0.9992996;0.8841292;0.4672424;0.03742076;0.9992996;0.26E-19]';
    f_opt = -0.8660254038;
elseif P== 79
    Name='HS 109';
    % 2 linear inequalities
    % 2 nonlinear inequalities
    % 6 nonlinear equalities
    b_L=[-0.55 -0.55]';
    b_U=inf*ones(2,1);
    A=[0 0 -1 1 0 0 0 0 0;
        0 0 1 -1 0 0 0 0 0]; 
    c_L = [0;0;0;0;0;0;0;0];
    c_U = [inf;inf;0;0;0;0;0;0];
    x_0 = [0;0;0;0;0;0;0;0;0];%% Not feasible
    x_L = [0;0;-0.55;-0.55;196;196;196;-400;-400];
    x_U = [inf;inf;0.55;0.55;252;252;252;800;800]; 
    x_min = [400;900;-1;-1;100;100;100;200;200];
    x_max = [800;1500;2;2;300;300;300;600;400];
    % Something wrong here!!!  x_opt not feasible
%     x_opt = [674.8881;1134.170;0.1335691;0.3711526;252.0000;252.0000;201.456;426.661;368.494]';
%     f_opt = 5326.06928;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [ 628.864531144928 1122.894687439662 0.096473078785 ...
               -0.219983362605 300 300 276.196437217777 ...
              270.332147445771 187.743505072573];
    f_opt = 5120.258595054240;
elseif P== 80
    Name='HS 111'; 
    b_L=[];b_U=[];A=[]; 
    c_L = [2;1;1];
    c_U = [2;1;1];
    x_0 = [-2.3;-2.3;-2.3;-2.3;-2.3;-2.3;-2.3;-2.3;-2.3;-2.3];%% Not feasible
    x_L = -100*ones(10,1);
    x_U = 100*ones(10,1); 
    x_min = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2];
    x_max = [2;2;2;2;2;2;2;2;2;2];
    x_opt = [-3.201212;-1.912060;-.2444413;-6.537489;-.7231524;-7.267738;-3.596711;-4.017769;-3.287462;-2.335582]';
    f_opt = -47.76109026;
elseif P== 81
    % 3 linear equalities
    Name='HS 112';
    b_L=[2;1;1];
    b_U=[2;1;1];
    A = [1 2 2 0 0 1 0 0 0 1;
        0 0 0 1 2 1 1 0 0 0;
        0 0 1 0 0 0 1 1 2 1];
    c_L = [];   
    c_U = [];
    x_0 = 0.1*ones(10,1); %Not feasible
    x_L = 1e-6*ones(10,1);
    x_U = []; 
    x_min = -ones(10,1);
    x_max = ones(10,1);
%     x_opt = [0.01773548;0.08200180;0.8825646;0.7233256E-3;0.4907851;0.4335469E-3;0.01727298;
%              0.007765639;0.01984929;0.05269826]';
%     f_opt = -47.707579;
    % Used glcCluster to improve f_opt.
    x_opt = [0.040311430138947 0.147808521731234 0.783262713386534 ...
             0.001414187429183 0.485252747699791 0.000693422931647 ...
             0.027386894239589 0.017937821913509 0.037279946883249 ...
             0.096852676693870];
    f_opt = -47.760726453844370;
elseif P== 82
    Name='HS 113';
    % 3 linear inequalities  
    % 5 nonlinear inequalities
    b_L=[-105 0 -12];
    b_U=inf*ones(3,1);
    A=[-4 -5 0 0 0 0 3 -9 0 0;
        -10 8 0 0 0 0 17 -2 0 0;
        8 -2 0 0 0 0 0 0 -5 2]; 
    c_L = [0;0;0;0;0];
    c_U = [];
    x_0 = [2;3;5;5;1;2;7;3;6;10];% Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10];
    x_opt = [2.171996;2.363683;8.773926;5.095984;0.9906548;1.430574;1.321644;9.828726;8.280092;8.37527]';
    f_opt = 24.30622091;
elseif P== 83
    Name='HS 114';
    % 4 linear inequalities
    % 1 linear equality
    % 3 nonlinear inequalities
    % 2 nonlinear equalities
    a=0.99;
    b=0.9;
    b_L=[-35.82 133 35.82 -133 0];
    b_U=[inf inf inf inf 0];
    A=[0 0 0 0 0 0  0  0  -b  -0.222;
        0 0 0 0 0 0  3  0   0  -0.99;
        0 0 0 0 0 0  0  0 (1/b) 0.222;
        0 0 0 0 0 0 -3  0   0   a;
        -1 0 0 1.22 -1  0   0   0 0 0]; 
    c_L = [0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;0;0];
    x_0 = [1745;12000;110;3048;1974;89.2;92.8;8;3.6;145];% Not feasible
    x_L = [0.00001;0.00001;0.00001;0.00001;0.00001;85;90;3;1.2;145];
    x_U = [2000;16000;120;5000;2000;93;95;12;4;162]; 
    x_min = [1000;10000;10;100;100;100;10;1;1;100];
    x_max = [2000;16000;70;3500;300;3000;200;200;3;200];
    %x_opt = [1698.096;15818.73;54.10228;3031.226;2000;90.11537;95;10.49336;1.561636;153.53535]';
    % New x_opt with slightly better f_opt
    x_opt = [1698.09473109545 15818.61065229285 54.10269773086 3031.22518942250 ...
             2000 90.11542432970 95 10.49329600172 1.56163636364 153.53535353535];
    f_opt = -1768.806963833686;
elseif P==84
    Name='HS 116';
    % 4  linear inequalities
    %10 nonlinear inequalities
    b_L=[0 -1 0 50];
    b_U=[inf inf inf 250];
    A=[0 -1 1 0 0 0 0 0 0 0 0 0 0;
        0  0 0 0 0 0 -0.002 0.002 0 0 0 0 0;
        -1 1 0 0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0 1 1 1]; 
    c_L = [0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [0.5;0.8;0.9;0.1;0.14;0.5;.489;.80;.650;.450;.150;.150;.150];% Not feasible
    x_L = [0.1;0.1;0.1;0.0001;0.1;0.1;0.1;0.1;500;0.1;1;0.0001;0.0001];
    x_U = [1;1;1;0.1;0.9;0.9;1000;1000;1000;500;150;150;150]; 
    x_min = [-1;-1;-1;-1;-1;-1;400;50;400;-1;10;50;-1];
    x_max = [2;2;2;2;2;2;600;100;600;2;21;80;2];
    % Something wrong here!!!  x_opt not feasible
%     x_opt = [0.8037703;0.8999860;0.9709724;0.9999952;0.1908154;0.4605717;574.0803;74.08043;500.0162;0.1;20.23413;77.34755;0.00673039]';
%     f_opt = 97.588409;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [0.8037731573767    0.8999858005409   0.9709830961312  0.1 ...
             0.1908131759963    0.4606603781675 574.0775733774242      ...
             74.0775733774241 500.0161605028427   0.1000000000000      ...
             20.2330906412816  77.3476900603028   0.0067288520660 ];
    f_opt = 97.587509553650435;
elseif P==85
    % 5 nonlinear inequalities
    Name='HS 117';
    b_L=[];b_U=[];A=[];
    c_L = [0;0;0;0;0];   
    c_U = [];
    x_0 = 0.001*[1 1 1 1 1 1 60000 1 1 1 1 1 1 1 1];%Feasible
    x_L = zeros(15,1);
    x_U = [];
    x_min =15*ones(15,1);
    x_max =12*ones(15,1);
%     x_opt = [0 0 5.174136 0 3.061093  11.83968 0 0 0.1039071 0 ...
%              0.2999929 0.3334709 0.3999910 0.4283145 0.2239607]';
%     f_opt = 32.348679;
    % New x_opt with slightly better f_opt
    x_opt = [0 0 5.174040715726261 0 3.061108657770277 11.839545808881381 ... 
             0 0 0.103896225611069 0 0.300000010139264  0.333467615710016 ...
                 0.399999993639061   0.428310101707230  0.223964867555319];
    f_opt = 32.348678965723366;
elseif P==86
    % 8  linear equalities
    Name='HS 119';
    b_L=[2.5  1.1 -3.1  -3.5  1.3  2.1  2.3  -1.5];
    b_U=[2.5  1.1 -3.1  -3.5  1.3  2.1  2.3  -1.5];
    A = [0.22  0.20  0.19  0.25  0.15  0.11  0.12  0.13   1.0  0.0  0.0 0.0 0.0 0.0 0.0 0.0;
        -1.46  0.0  -1.30  1.82 -1.15  0.0   0.80  0.0    0.0  1.0 0.0 0.0 0.0 0.0 0.0 0.0;
        1.29 -0.89  0.0   0.0  -1.16 -0.96  0.0  -0.49   0.0  0.0  1.0 0.0 0.0 0.0 0.0 0.0;
        -1.10 -1.06  0.95 -0.54  0.0  -1.78 -0.41  0.0    0.0  0.0  0.0 1.0 0.0 0.0 0.0 0.0;
        0.0   0.0   0.0  -1.43  1.51  0.59 -0.33 -0.43   0.0  0.0  0.0 0.0 1.0 0.0 0.0 0.0; 
        0.0  -1.72 -0.33  0.0   1.62  1.24  0.21 -0.26   0.0  0.0  0.0 0.0 0.0 1.0 0.0 0.0;
        1.12  0.0   0.0   0.31  0.0   0.0   1.12  0.0   -0.36 0.0  0.0 0.0 0.0 0.0 1.0 0.0;
        0.0   0.45  0.26 -1.10  0.58  0.0  -1.03  0.10   0.0  0.0  0.0 0.0 0.0 0.0 0.0 1.0];          
    c_L = [];   
    c_U = [];
    x_0 = 10*ones(16,1); %Not feasible
    x_L = zeros(16,1);
    x_U = 5*ones(16,1); 
    x_min = -1*ones(16,1);
    x_max = 10*ones(16,1);
    % Something wrong here!!!  x_opt not feasible, chs_f(x_opt) ~= f_opt
%     x_opt = [0.3984735;0.7919832;0.2028703;0.8443579;1.126991;0.9347387;...
%             1.681962;1.553009;1.5667870;0;0;0;0.6602041;0;0.6742559;0]';
%     f_opt = 244.899698;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [0.039848887485877 0.791984373292997 0.202868317085505 ...
             0.844360022151986 1.269909934712237 0.934735378604400 ...
             1.681961715678740 0.155300989591008 1.567869467928664 ...
             0 0 0 0.660203748583386 0 0.674253526042831 0];
    f_opt = 244.8996975176148;
elseif P==87
    Name='HS 215';
    %1 nonlinear inequality
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1;1]; %Feasible
    x_L = [0;0];
    x_U = []; 
    x_min = [-2;-2;];
    x_max = [2;2];
    x_opt = [0;0]';
    f_opt = 0;
elseif P==88
    Name='HS 216';
    %1 nonlinear equality
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [-1.2;1]; %Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-10;-10;];
    x_max = [10;10];
    x_opt = [2;4]';
    f_opt = 1;
elseif P==89
    Name='HS 217';
    %1 linear inequality
    %1 nonlinear equality
    b_L=-1;b_U=inf;A=[1 -2]; 
    c_L = 0;
    c_U = 0;
    x_0 = [10;10]; %Not feasible
    x_L = [0;-inf];
    x_U = [inf;inf]; 
    x_min = [-2;-2];
    x_max = [2;2];
    x_opt = [0.6;0.8]';
    f_opt = -0.8;
elseif P==90
    Name='HS 218';
    %1 nonlinear equality
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [9;100]; %Feasible
    x_L = [-inf;0];
    x_U = []; 
    x_min = [-2;-2];
    x_max = [2;2];
    x_opt = [0;0]';
    % This is not correct. f(x_opt) = 0.
    %f_opt = 100;
    f_opt = 0;
end

c  = 'chs_c';
dc = 'chs_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('chs_f','chs_g','chs_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P      = P;

if P==49 
    Prob.FUNCS.dc = '';
end

function Prob = chs_prob3(P)

HessPattern = []; pSepFunc = []; ConsPattern = []; f_Low = [];

if P==91
    Name='HS 219';
    %2 nonlinear inequalities
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [10;10;10;10]; %Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-2;-2;-2;-2];
    x_max = [2;2;2;2];
    x_opt = [1;1;0;0]';
    f_opt = -1;
elseif P==92
    Name='HS 220';
    % 1 nonlinear equality   
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2500;2500]; %Feasible
    x_L = [1;0];
    x_U = []; 
    x_min = [-20000;-20000];
    x_max = [30000;30000];
    x_opt = [1;0]';
    f_opt = 1;
elseif P==93
    Name='HS 221';
    % 1 nonlinear inequality   
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [0.25;0.25]; %Feasible
    x_L = [0;0];
    x_U = []; 
    x_min = [-2;-2];
    x_max = [2;2];
    x_opt = [1;0]';
    f_opt = -1;
elseif P==94
    Name='HS 222';
    % 1 nonlinear inequality   
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1.3;0.2]; % Not feasible
    x_L = [0;0];
    x_U = []; 
    x_min = [-3;-3];
    x_max = [3;3];
    x_opt = [1.5;0]';
    % This is not correct. f(x_opt) = -1.5
    f_opt = -1.5;
elseif P==95
    Name='HS 223';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0;0]; % Feasible
    x_L = [0;0];
    x_U = [10;10]; 
    x_min = [-15;-15];
    x_max = [15;15];
    x_opt = [0.834;10]';
    f_opt = -0.834032;
elseif P==96
    % Many strange things here!!!
    % Should it be equalities?? Then infeasible..
    % If constraints correct, then wrong optimal solution.
    % Maybe x_L = [0,0], then optimum correct??
    % Otherwise infinately many solutions:
    % {x: x1+x2=1, x1 < -(1+sqrt(5))/2 or x2 < -(1+sqrt(5))/2 }
    % 5 nonlinear equalities
    Name='HS 225';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0];
    c_U = [];
    x_0 = [3;1]; % Feasible
    x_L = [];
    x_U = []; 
    x_min = [-15;-15];
    x_max = [15;15];
    x_opt = [1;1]';
    f_opt = 2;
elseif P==97
    % One inequality redundant!!
    % 2 nonlinear inequalities
    Name='HS 226';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0.8;0.05]; % Feasible
    x_L = [0;0];
    x_U = []; 
    x_min = [-15;-15];
    x_max = [15;15];
    % Something wrong here!!!  x_opt not feasible, chs_f(x_opt) ~= f_opt
%     x_opt = [1;1]';
%     f_opt = 2;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [1/sqrt(2) 1/sqrt(2)];
    f_opt = -0.5;
elseif P==98
    % 2 nonlinear inequalities
    Name='HS 227';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0.5;0.5]; % Feasible
    x_L = [];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [1;1]';
    f_opt = 1;
elseif P==99
    % 2 nonlinear inequalities
    Name='HS 228';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0;0]; % Feasible
    x_L = [];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [0;-3]';
    f_opt = -3;
elseif P==100
    % 2 nonlinear inequalities
    Name='HS 230';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [0;0]; %Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [0.5;0.375]';
    f_opt = 0.375;
elseif P==101
    % 2 nonlinear inequalities
    Name='HS 231';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [];
    x_0 = [-1.2;1]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [1;1]';
    f_opt = 0;
elseif P==102
    % 3 linear inequalities
    Name='HS 232';
    b_L=[0 0 -6]';
    b_U=inf*ones(3,1);
    A=[1/sqrt(3)  -1;1  sqrt(3);-1  -sqrt(3)]; 
    c_L = [];
    c_U = [];
    x_0 = [2;0.5]; %Feasible
    x_L = [0;0];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [3;1.732]';
    f_opt = -1;
elseif P==103
    % 1 nonlinear inequality
    Name='HS 233';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [1.2;1]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [1;1]';
    f_opt = 0;
elseif P==104
    % 1 nonlinear inequality
    Name='HS 234';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = [];
    x_0 = [0;0]; %Not feasible
    x_L = [0.2;0.2];
    x_U = [2;2]; 
    x_min = [-5;-5];
    x_max = [5;5];
    x_opt = [0.2;0.2]';
    f_opt = -0.8;
elseif P==105
    % 1 nonlinear equality
    Name='HS 235';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [-2;3;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1];
    x_max = [1;1;1];
    x_opt = [-1;1;0]';
%     f_opt = -0.04;   % Not possible, f = sum of squares
    f_opt = 0.04;
elseif P==106
    Name='HS 236';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [90;10];
    x_L = [0;0]; 
    x_U = [75;65];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [75;65]';
%     f_opt = -58.9034;
    f_opt = -58.903436009034657;
elseif P==107
    Name='HS 237';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [95;10];
    x_L = [54;-inf]; 
    x_U = [75;65];
    %x_min = [-2 -2];
    x_min = [0 0];
    x_max = [2 2];
    x_opt = [75;65]';
%     f_opt = -58.9034;
    f_opt = -58.903436009034657;
elseif P==108
    Name='HS 238';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [95;10];
    x_L = [-inf;-inf]; 
    x_U = [75;65];
    %x_min = [-2 -2];
    x_min = [0 0];
    x_max = [2 2];
    x_opt = [75;65]';
%     f_opt = -58.9034;
    f_opt = -58.903436009034657;
elseif P==109
    Name='HS 239';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [95;10];
    x_L = [0;0]; 
    x_U = [75;65];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [75;65]';
%     f_opt = -58.9034;
    f_opt = -58.903436009034657;
elseif P==110
    Name='HS 248';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;0];
    x_0 = [-.1;-1;.1];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [.6;.8;0]';
    f_opt = -.8;
elseif P==111
    Name='HS 249';
    b_L=[];b_U=[];A=[];
    c_L = 0;
    c_U = [];
    x_0 = [1;1;1];
    x_L = [1;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Infiniately many solutions!! {x: x1^2 + x2^2 = 1, x3=0 }
    x_opt = [1;0;0]';
    f_opt = 1;
elseif P==112
    Name='HS 250';
    b_L=[0;-72];b_U=[inf;inf];A=[1 2 2;-1 -2 -2]; 
    c_L = [];
    c_U = [];
    x_0 = [10;10;10];
    x_L = [0;0;0]; 
    x_U = [20;11;42];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [20;11;15]';
    f_opt = -3300;
elseif P==113
    Name='HS 251';
    b_L=-72;b_U=inf;A=[-1 -2 -2]; 
    c_L = [];
    c_U = [];
    x_0 = [10;10;10];
    x_L = [0;0;0]; 
    x_U = [42;42;42];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Something wrong here!!!  chs_f(x_opt) ~= f_opt
    %x_opt = [20;12;12]';
    x_opt = [24;12;12]';
    f_opt = -3456;
elseif P==114
    Name='HS 252';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [-1;2;2];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [-1;1;0]';
    f_opt = .04;
elseif P==115
    Name='HS 253';
    b_L=-inf;b_U=30;A=[3 0 3]; 
    c_L = [];
    c_U = [];
    x_0 = [0;2;0];
    x_L = [0;0;0]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Something strange here!! If all coefficients are correct,
    % then there exists better solutions.
%     x_opt = [.3333;.3333;.3333]';
%     f_opt = 87.3794;
    x_opt = [5 5 5];
    f_opt = 69.282032302755098;
elseif P==116
    Name='HS 254';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [1;1;1];
    x_L = [-inf;-inf;1]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [0;1.732;1]';
    f_opt = -1.73205;
elseif P==117
    Name='HS 263';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0];
    c_U = [inf;inf;0;0];
    x_0 = [10;10;10;10];
    x_L = [-inf;-inf;-inf;-inf]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    x_opt = [1;1;0;0]';
    f_opt = -1;
elseif P==118
    Name='HS 264';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [0;0;0;0];
    x_L = [-inf;-inf;-inf;-inf]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    x_opt = [0;1;2;-1]';
    f_opt = -44;
elseif P==119
    Name='HS 265';
    b_L=[1;1];b_U=[1;1];A=[1 1 0 0;0 0 1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [0;0;0;0];
    x_L = [0;0;0;0]; 
    x_U = [1;1;1;1];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    x_opt = [1;0;1;0]';
    f_opt = .974747;
elseif P==120
    Name='HS 270';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [1.1;2.1;3.1;4.1;-1];
    x_L = [1;2;3;4;-inf]; 
    x_U = [inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2];
    x_opt = [1 2 3 4 2];
    f_opt = -1;
elseif P==121
    Name='HS 284';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
    x_opt = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]';
    f_opt = -1840;
elseif P==122
    Name='HS 285';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    x_L = []; 
    x_U = [];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2 2 2 2 2 2 2];
    x_opt = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]';
    f_opt = -8252;
elseif P==123
    Name='HS 315';
    b_L=-1;b_U=inf;A=[1 -2]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [-0.1;-0.9];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [.6;.8]';
    f_opt = -.8;
elseif P==124
    Name='HS 316';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-10 -10];
    x_max = [10 10];
    x_opt = [7.071;-7.071]';
    f_opt = 334.315;
elseif P==125
    Name='HS 317';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-10 -10];
    x_max = [10 10];
    x_opt = [7.352;-5.423]';
    f_opt = 372.467;
elseif P==126
    Name='HS 318';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [0 -10];
    x_max = [20 10];
    x_opt = [7.809;-3.748]';
    f_opt = 412.750;
elseif P==127
    Name='HS 319';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-5 -5];
    x_max = [15 15];
    x_opt = [8.492;-2.112]';
    f_opt = 452.404;
elseif P==128
    Name='HS 320';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-5 -5];
    x_max = [15 15];
    x_opt = [9.395;-.6846]';
    f_opt = 485.531;
elseif P==129
    Name='HS 321';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-5 -5];
    x_max = [15 15];
    x_opt = [9.816;-.1909]';
    f_opt = 496.112;
elseif P==130
    Name='HS 322';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-5 -5];
    x_max = [15 15];
    x_opt = [9.998;-.1999e-2]';
    f_opt = 499.96;
elseif P==131
    Name='HS 323';
    b_L=-2;b_U=inf;A=[1 -1]; 
    c_L = 0;
    c_U = inf;
    x_0 = [0;0];
    x_L = [0;0]; 
    x_U = [inf;inf];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [.5536;1.306]';
    f_opt = 3.79894;
elseif P==132
    Name='HS 324';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [2;2];
    x_L = [0;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-20 -20];
    x_max = [20 20];
    % Alternative optima!
    x_opt = [ 15.811388114194754  1.581138848747321 ;
             -15.811388114194754 -1.581138848747321];
    f_opt = 5;
elseif P==133
    Name='HS 325';
    b_L=-inf;b_U=1;A=[1 1]; 
    c_L = [0;0];
    c_U = [inf;0];
    x_0 = [-3;0];
    x_L = [-inf;-inf]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [-5 -5];
    x_max = [5 5];
    % New x_opt that satisfies constraint better. slightly worse f_opt 
%     x_opt = [-2.372;-1.836]';
%     f_opt = 3.79134;
    x_opt = [-2.372281335897289  -1.836377311044883];
    f_opt = 3.791341425601744;
elseif P==134
    Name='HS 326';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [4;3];
    x_L = [0;0]; 
    x_U = [inf;inf];
%     x_min = [-2 -2];
%     x_max = [2 2];
    x_min = [0 0];
    x_max = [10 10];
    % New x_opt that satisfies constraint better. slightly worse f_opt 
%     x_opt = [5.240;3.746]';
%     f_opt = -79.8078;
    x_opt = [5.239609115511264 3.746037752432500];
    f_opt = -79.807820846506957;
elseif P==135
    Name='HS 327';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [.42;5];
    x_L = [0.4;0.4]; 
    x_U = [inf;inf];
    x_min = [-2 -2];
    x_max = [5 5];
    x_opt = [.4219;3.746]';
    % This is not correct. f(x_opt) = 0.03063
    %f_opt = -79.8078;
    f_opt = 0.030631492380913;
end

% Define the Prob
c  = 'chs_c';
dc = 'chs_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('chs_f','chs_g','chs_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P      = P;

function Prob = chs_prob4(P)

HessPattern = []; pSepFunc = []; ConsPattern = []; f_Low = [];

if P==136
    Name='HS 329';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [14.35;8.6];
    % Something is wrong. x_L(1) = 13 but still x_opt(1) < 13???
    %x_L = [13;0]; 
    x_L = [10;0]; 
    x_U = [16;15];
    x_min = [-2 -2];
    x_max = [2 2];
    % Something is wrong here!!  chs_c(x_opt) not feasible. 
    % Small typo in x_opt(1). Change from 12.09 to 14.095
%     x_opt = [14.095;.8430]';
%     f_opt = -6961.81;
    % Use solution from glcDirect:
    x_opt = [14.095 0.842961];
    f_opt = -6961.813643511347;
elseif P==137
    Name='HS 330';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [2.5;2.5];
    x_L = [0;0]; 
    x_U = [5;5];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [1.287;.5305]';
    f_opt = 1.62058;
elseif P==138
    Name='HS 331';
    b_L=-inf;b_U=1;A=[1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [.5;.1];
    x_L = [.0001;.0001]; 
    x_U = [inf;1];
    x_min = [-2 -2];
    x_max = [2 2];
    x_opt = [.6175;.1039]';
    f_opt = 4.258;
elseif P==139
    % THIS PROBLEM HAS NO FEASIBLE SOLUTION WITHIN x_L and x_U!!!
    Name='HS 332';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [.75;.75];
    x_L = [0;0]; 
    x_U = [1.5;1.5];
    x_min = [-2 -2];
    x_max = [2 2];
    % This is not correct. f(x_opt) = 29.92
    % Also, x_opt is not feasible!!
    %f_opt = 114.95;
    x_opt = [.9114;.02928]';
    f_opt = 29.924379392278784;
elseif P==140
    Name='HS 335';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [0;0];
    x_0 = [1;1;1];
    x_L = [-inf;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [.2031e-5;.4472e-2;.2000e-2]';
    f_opt = -.447214e-2;
elseif P==141
    Name='HS 336';
    b_L=6;b_U=6;A=[5 5 -3]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0;0];
    x_L = [-inf;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [.5346;.534;-.2191]';
    f_opt = -.337896;
elseif P==142
    Name='HS 337';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [1;1;1];
    x_L = [-inf;1;-inf]; 
    x_U = [inf;inf;1];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [.5774;1.732;-.2026e-5]';
    f_opt = 6;
elseif P==143
    Name='HS 338';
    b_L=1;b_U=1;A=[.5 1 1]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0;0];
    x_L = [-inf;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 3 2];
    % Found better solution with glcCluster.
%     x_opt = [.3669;2.244;-1.427]';
%     f_opt = -7.20570;
    x_opt = [1.157072134057926 -1.578536067028963 2];
    f_opt = -7.830592038324628;
elseif P==144
    Name='HS 339';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [1;1;1];
    x_L = [0;0;0]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [3 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [2.380;.3162;1.943]';
%     f_opt = 3.36168;
    x_opt = [2.380 0.31568 1.943];
    f_opt = 3.361680361960383;
elseif P==145
    % Something is wrong here!!!
    % This is a local optimum, found by snopt.
    % But as x_L = x_U = [], an unbounded global solution exists.
    Name='HS 340';
    b_L=-inf;b_U=1.8;A=[1 2 2]; 
    c_L = [];
    c_U = [];
    x_0 = [1;1;1];
    x_L = [-inf;-inf;-inf]; 
    x_U = [1;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [.6;.3;.3]';
    f_opt = -.054;            
elseif P==146
    Name='HS 341';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [1;1;1];
    x_L = [0;0;0]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [4 4 4];
    x_opt = [4;2.828;2]';
    f_opt = -22.6274;
elseif P==147
    Name='HS 342';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = inf;
    x_0 = [100;100;100];
    x_L = [0;0;0]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [4 4 4];
    x_opt = [4;2.828;2]';
    f_opt = -22.6274;
elseif P==148
    Name='HS 343';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [22.3;.5;125];
    x_L = [0;0;0]; 
    x_U = [36;5;125];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [16.51;2.477;124]';
%     f_opt = -5.68478;
    % Also found alternative optima, which actually is tiny tiny better!
    x_opt = [16.466180355867100 2.489533668773100 124.312312115096100 ;
             26.614576872000207 0.952937989749282  76.910820671366878];
    f_opt = -5.684782500034809;
elseif P==149
    Name='HS 344';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2;2];
    x_L = [-inf;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.105;1.197;1.535]';
%     f_opt = .0325682;
    x_opt = [1.104859022986986 1.196674169558591 1.535262262103906];
    f_opt = 0.032568200255028;
elseif P==150
    % Identical to P=149!!! ??? (also for chs_f and chs_c)
    Name='HS 345';
    b_L=[];b_U=[];A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [2;2;2];
    x_L = [-inf;-inf;-inf]; 
    x_U = [inf;inf;inf];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.105;1.197;1.535]';
%     f_opt = .0325682;
    x_opt = [1.104859022986986 1.196674169558591 1.535262262103906];
    f_opt = 0.032568200255028;
elseif P==151
    Name='HS 346';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [22.3;.5;125];
    x_L = [0;0;0]; 
    x_U = [36;5;125];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [16.51;2.477;124]';
%     f_opt = -5.68478;
    % Also found alternative optima, which actually is tiny tiny better!
    x_opt = [16.4661803558671 2.4895336687731 124.3123121150961 ;
             18.7761992001614 1.9146456288525 109.0182803890518];
    f_opt = -5.684782500032257;
elseif P==152
    Name='HS 347';
    b_L=1;b_U=1;A=[1 1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [.7;.2;.1];
    x_L = [0;0;0]; 
    x_U = [1;1;1];
    x_min = [-2 -2 -2];
    x_max = [2 2 2];
    x_opt = [0;0;1]';
    f_opt = 17374.6;
elseif P==153
    % Something is wrong here!!!
    % If P.b_U is correct, then x_opt is infeasible.
    % If P.b_U(2) is changed from 0 to 1, then everything is ok.
    Name='HS 353';
    b_L=[5;0];
%     b_U=[inf;0];
    b_U=[inf;1];
    A=[2.3 5.6 11.1 1.3;
        1 1 1 1]; 
    c_L = 0;
    c_U = inf;
    x_0 = [0;0;.4;.6];
    x_L = [0;0;0;0]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    x_opt = [0;0;.3776;.6224]';
    f_opt = -39.9337;
elseif P==154
    Name='HS 354';
    b_L=1;
    b_U=inf;
    A=[1 1 1 1]; 
    c_L = [];
    c_U = [];
    x_0 = [3;-1;0;1];
    x_L = [-inf;-inf;-inf;-inf]; 
    x_U = [20;20;20;20];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    % Probably typo!!
    %x_opt = [.5034;-.4557;.2358;.3064]';
    x_opt = [.5034;-.04557;.2358;.3064]';
    f_opt = .113784;
elseif P==155
    Name='HS 355';
    b_L=[];b_U=[]; A=[]; 
    c_L = 0;
    c_U = 0;
    x_0 = [0;0;0;0];
    x_L = [.1;.1;0;0]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.917;.1;0;1.972]';
%     f_opt = 69.6755;
    % Also found alternative optima, which actually is tiny tiny better!
    x_opt = [1.916633142925654 0.1 0 1.971811696396848];
    f_opt = 69.675463076595392;
elseif P==156
    Name='HS 356';
    b_L=0;b_U=inf; A=[-1 0 0 1]; 
    c_L = [0;0;0;0];
    c_U = [inf;inf;inf;inf];
    x_0 = [1;7;8;1];
    x_L = [.125;0;0;-inf]; 
    x_U = [inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2];
    x_max = [2 2 2 2];
    % Found better solution with glcCluster.
%     x_opt = [.2444;6.219;8.291;.2444]';
%     f_opt = 2.38116;
    x_opt = [0.154541563495586 16.368103811468124 ...
             3.007055191781054 0.322931111645804];
    f_opt = 1.850599229908934;
elseif P==157
    Name='HS 360';
    b_L=[];b_U=[]; A=[]; 
    c_L = [0;0];
    c_U = [inf;inf];
    x_0 = [2.52;2;37.5;9.25;6.8];
    x_L = [0;1.2;20;9;6.5]; 
    x_U = [inf;2.4;60;9.3;7];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [5 5 60 10 10];
    x_opt = [4.537;2.4;60;9.3;7]';
    f_opt = -.528034e7;
elseif P==158
    Name='HS 361';
    b_L=[];b_U=[]; A=[]; 
    c_L = [0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf];
    x_0 = [2.52;2;37.5;9.25;6.8];
    x_L = [0;1.2;20;9;-inf]; 
    x_U = [inf;2.4;60;9.3;7];
    x_min = [-2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2];
    x_opt = [.6813;.24;20;9.3;7]';
    f_opt = -776412.12;
elseif P==159
    Name='HS 365';
    b_L=[];b_U=[]; A=[]; 
    c_L = [0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf];
    x_0 = [3;0;2;-1.5;1.5;5;0];
    x_L = [0;-inf;0;-inf;1;-inf;1]; 
    x_U = [inf;inf;inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2];
    x_opt = [4.83;0;4.83;1;2.41;2.41;1]';
    f_opt = 23.3137;
elseif P==160
    Name='HS 366';
    b_L = [-inf;-inf]; b_U = [1;1]; 
    A   = [     0 0      0 0 .022556 0 -.007595 ;
           -.0005 0 .00061 0 0       0        0]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [1745;110;3048;89;92.8;8;145];
    x_L = [1;1;1;85;90;3;145]; 
    x_U = [2000;120;5000;93;95;12;162];
    x_min = [-2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2];
    % x_opt not feasible!!
%     x_opt = [905.4;36.39;2381.48;88.99;95;12;153.54]';
%     f_opt = 704.306;
    % Use glcCluster to find feasible solution
    x_opt = [1698.185162338148 53.664931346193 3031.299313391925 ...
               90.109752694742 95 10.499428077951 153.535354545455 ];
    f_opt = 1227.225195121009;
elseif P==161
    Name='HS 367';
    b_L = [-inf;-inf;5]; 
    b_U = [10;5;5];
    A   = [ 1 1 1 1 1  1 1 ;
            1 1 1 1 0  0 0 ;
            0 0 0 2 1 .8 1]; 
    c_L = [0;0];
    c_U = [inf;0];
    x_0 = [.1;.1;.1;.1;.1;.1;.1];
    x_L = [0;0;0;0;0;0;0]; 
    x_U = [inf;inf;inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.47;1.98;.35;1.2;.57;.78;1.41]';
%     f_opt = -37.4130;
    x_opt = [1.468809853711926 1.983971313449347 0.351877282326562 ...
             1.195341550512166 0.569398237110642 0.784745725176464 ...
             1.412122081723855];
    f_opt = -37.412959541261586;
elseif P==162
    Name='HS 369';
    b_L=[-inf;-inf;-inf];b_U=[1;1;1]; 
    A=[0 0 0 .0025 0 .0025 0 0;
        0 0 0 -.0025 .0025 0 .0025 0;
        0 0 0 0 -.01 0 0 .01]; 
    c_L = [0;0;0];
    c_U = [inf;inf;inf];
    x_0 = [5000;5000;5000;200;350;150;225;425];
    x_L = [100;1000;1000;10;10;10;10;10]; 
    x_U = [1e4;1e4;1e4;1e3;1e3;1e3;1e3;1e3];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2];
    x_opt = [579.3;1360;5101;182;295.6;217.9;286.4;395.6]';
    f_opt = 7049.24;
elseif P==163
    Name='HS 372';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [300;-100;-.1997;-127;-151;379;421;460;426];
    x_L = [-inf;-inf;-inf;0;0;0;0;0;0]; 
    x_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2];
    % Change x_opt somewhat to improve feasability
%     x_opt = [523.31;-156.95;-.2;29.61;86.62;47.33;26.24;22.92;39.47]';
%     f_opt = 13390.1;
    x_opt = [523.3055381965517 -156.9478431234874  -0.1996645694891 ...
              29.6080104373793   86.6155546684712  47.3267002003196 ...
              26.2355969860397   22.9159847176167  39.4707368850624];
    f_opt = 13390.09311947957;
elseif P==164
    Name='HS 373';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;];
    c_U = [0;0;0;0;0;0;];
    x_0 = [300;-100;-.1997;-127;-151;379;421;460;426];
    x_L = []; 
    x_U = [];
    x_min = [-2   -200 -2  -2 -100 -2 -2 -2 -40];
    x_max = [ 600    2  2  30    2 50 30 30   2];
    % Something wrong here!!!  x_opt not feasible and chs_f(x_opt) ~= f_opt
%     x_opt = [523.31;-156.95;-.2;29.61;-86.62;47.33;26.24;22.92;-39.47]';
%     f_opt = 13390.1;
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [523.3055409521729 -156.9478460145828  -0.1996645670108 ...
              29.6080102497378  -86.6155542855247  47.3267004997901 ...
              26.2355969167998   22.9159841912966 -39.4707378586784];
    f_opt = 13390.09311947915;
elseif P==165
    Name='HS 374';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;
        inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [.1;.1;.1;.1;.1;.1;.1;.1;.1;.1];
    x_L = [-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf]; 
    x_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2 2];
    % Found alternative optima!
    x_opt = [ 0.2182  0.2326  0.2785  0.2681  0.212  0.1259  0.0341 -0.0261 -0.1422 0.2333 ;
             -0.2182 -0.2326 -0.2785 -0.2681 -0.212 -0.1259 -0.0341  0.0261  0.1422 0.2333 ];
    f_opt = .233264;
elseif P==166
    Name='HS 375';
    A    = ones(8,10);
    for i=1:8, A(i,i) = 1/2; end     
    b_L   = ones(8,1);
    b_U   = ones(8,1);
    c_L   = 4;
    c_U   = 4;
    x_0   = ones(10,1);
    x_L   = [];
    x_U   = [];
    x_min = -ones(10,1);
    x_max =  ones(10,1);
    x_opt = [0.1064 0.1064 0.1064 0.1064 0.1064 0.1064 0.1064 0.1064 2.843 -2.642];
    f_opt = -15.16;
elseif P==167
    Name='HS 376';
    b_L=.255;b_U=.255;A=[0 0 0 0 0 0 0 0 1 1]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    c_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_0 = [10;.005;.0081;100;.0017;.0013;.0027;.002;.15;.105];
    x_L = [0;0;0.5e-4;10;.5e-4;.5e-4;.5e-4;.5e-4;.5e-4;.5e-4]; 
    x_U = [10;.1;.0081;1000;.0017;.0013;.0027;.002;1;1];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 650 2 2 2 2 2 2];
    x_opt = [.14727;.1;.0081;628.72;.0017;.001182;.0027;.00135;.15741;.097593]';
    f_opt = -4430.88;
    ConsPattern = spalloc(14,10,48);
    ConsPattern(1, [1,3,4])  = 1;
    ConsPattern(2, [1,4,5,9])  = 1;
    ConsPattern(3, [1,4,6,10])  = 1;
    ConsPattern(4, [1,4,7])  = 1;
    ConsPattern(5, [1,4,8])   = 1;
    ConsPattern(6, [2,4,5,9])  = 1;
    ConsPattern(7, [2,4,6,10])  = 1;
    ConsPattern(8, [2,4,7])  = 1;
    ConsPattern(9, [2,4,8])  = 1;
    ConsPattern(10,[2,3,4]) = 1;
    ConsPattern(11,[2,3,4]) = 1;
    ConsPattern(12,[2,3,4]) = 1;
    ConsPattern(13,[2,3,4]) = 1;
    ConsPattern(14,[3,4,5,6,8])  = 1;
    
elseif P==168
    Name='HS 377';
    b_L=[2;1;1];b_U=[2;1;1];
    A=[1 -2 2 0 0 1 0 0 0 1;
        0 0 0 1 -2 1 1 0 0 0;
        0 0 1 0 0 0 1 1 2 1]; 
    c_L = [];
    c_U = [];
    x_0 = [.1;.1;.1;.1;.1;.1;.1;.1;.1;.1];
    a = .1e-4;
    x_L = [a;a;a;a;a;a;a;a;a;a]; 
    x_U = [10;10;10;10;10;10;10;10;10;10];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2 2];
    x_opt = [10;10;1;10;9.5;10;1e-4;1e-4;1e-4;1e-4]';
    f_opt = -795.001;
elseif P==169
    Name='HS 378';
    b_L=[];b_U=[];A=[]; 
    c_L = [2;1;1];
    c_U = [2;1;1];
    x_0 = -2.3*ones(10,1);
    x_L = [-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf;-inf]; 
    x_U = [inf;inf;inf;inf;inf;inf;inf;inf;inf;inf];
    x_min = [-2 -2 -2 -2 -2 -2 -2 -2 -2 -2];
    x_max = [2 2 2 2 2 2 2 2 2 2];
    % Found alternative optima!
    x_opt = [-3.2024 -1.9123 -0.2444 -15.67 -0.7217 -7.2736 -3.5965 -4.0206 -3.2885 -2.3344 ;
             -3.2024 -1.9123 -0.2444 -10.00 -0.7217 -7.2736 -3.5965 -4.0206 -3.2885 -2.3344 ];
    f_opt = -47.7610;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif P==170
    % 3 nonlinear inequality
    Name='HS 380';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0];
    c_U = [];
    x_0 = [4;4;4;4;4;4;4;4;4;4;4;4]; %Not feasible
    x_L = [0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1];
    x_U = [100;100;100;100;100;100;100;100;100;100;100;100]; 
    x_min = [0;0;0;0;0;0;0;0;0;0;0;0];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10];
    % Something wrong here!!!  chs_f(x_opt) = 99805.76450630033
%     f_opt = 3.16859;
%    x_opt = [2.6632;4.5173;7.1338;2.2373;4.0784;1.3183;4.1252;2.8562;1.6756;2.1789;5.1234;6.6593]';
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [4 4 4 1 1 4 4 4 4 4 100 4];
    f_opt = 99453.48361219863;
elseif P==171
    % 3 nonlinear inequalities
    % 1 linear equality
    Name='HS 382';
    b_L=1;b_U=1;A=[1 1 1 1 1 1 1 1 1 1 1 1 1]; 
    c_L = [0;0;0];
    c_U = [];
    x_0 = [0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1;0.1]; %Feasible
    x_L = [0;0;0;0;0;0;0;0;0;0;0;0;0];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10];
    x_opt = [0.132;0;0;0;0;0.3263;0;0;0.5167;0;0;0;0.025]';
    f_opt = 1.03831;
elseif P==172
    % 1 linear equality
    Name='HS 383';
    b_L=1;b_U=1;
    A=[5.47934 0.83234 0.94749 1.11082 2.64824 1.55868 1.73215 3.90896 2.74284 2.60541 5.96184 3.29522 1.83517 2.81372]; 
    c_L = [];
    c_U = [];
    x_0 = [0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01;0.01];%Not Feasible
    x_L = [0;0;0;0;0;0;0;0;0;0;0;0;0;0];
    x_U = 10*[0.04;0.04;0.04;0.04;0.04;0.03;0.03;0.03;0.03;0.03;0.03;0.03;0.03;0.03]; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [1;1;1;1;1;1;1;1;1;1;1;1;1;1];
    % Something wrong here!!!  chs_f(x_opt) = 728592.2676748921, ~= P.f_opt
%     f_opt = 728565;
%     x_opt = [0.04;0.0382;0.033;0.0358;0.0303;0.0279;0.0265;0.0249;0.023;0.0216;0.0202;0.0192;0.0203;0.0253]';
    % Used glcCluster to find new feasible x_opt AND better f_opt.
    x_opt = [0.058372065541090 0.033283416185679 0.031195444831010  ...
             0.028810871415281 0.026383300539696 0.024314859157350  ...
             0.023062893246090 0.021707306750786 0.020070984161143  ...
             0.018797409243888 0.017571863328282 0.016709563269051  ...
             0.017677123424725 0.022047799775992];
    f_opt = 687865.6575787148;
elseif P==173
    % 10 nonlinear inequality
    Name='HS 384';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    % Change x_opt somewhat to improve feasability
%     x_opt = [0.861;0.917;0.92;0.896;1.037;0.973;0.822;1.2;1.156;1.144;1.03;0.909;1.082;0.847;1.172]';
%     f_opt = -8310.26;
    x_opt = [0.860953828013186 0.917361283651153 0.919736443049539 ...
             0.896005632656558 1.037294875168159 0.973089007736045 ...
             0.822436319111331 1.198721750529998 1.156334977108481 ...  
             1.144386760299319 1.030568064307347 0.909494608001779 ...
             1.082045210893987 0.846823767791531 1.172372176008728];
    f_opt = -8310.258974248296;
elseif P==174
    % 10 nonlinear inequality
    Name='HS 385';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    % x_opt not feasible!!
%     x_opt = [0.813;0.113;1.086;0.998;1.075;1.069;0.628;1.093;0.914;0.862;1.005;0.887;0.987;1.041;1.186]';
%     f_opt = -8315.29;
    % Use glcCluster to find feasible solution
    x_opt = [0.8135 1.1328 1.0861 0.9983 1.0755 1.0689 0.6278 1.0930 ...
             0.9136 0.8619 1.0047 0.8774 0.9867 1.0411 1.1861];
    f_opt = -8315.214;
elseif P==175
    % 11 nonlinear inequality
    Name='HS 386';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.004;1.087;1.103;1.031;0.929;1.257;0.761;0.857;1.09;0.981;0.851;0.966;0.9060;0.838;0.809]';
%     f_opt = -8154.37;
    x_opt = [1.0043 1.0871 1.1034 1.0307 0.9286 1.2568 0.7606 0.8569 ...
             1.0898 0.9812 0.8511 0.9656 0.9064 0.8380 0.8093];
    f_opt = -8164.3825;
elseif P==176
    % 11 nonlinear inequality
    Name='HS 387';
    b_L=[];b_U=[];A=[]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    % Change x_opt somewhat to improve feasability
%     x_opt = [1.013;1.016;1.031;0.997;0.985;1.037;0.993;0.972;0.999;0.995;0.969;1.008;0.982;0.991;0.978]';
%     f_opt = -8250.14;
    x_opt = [1.0125 1.0159 1.0309 0.9970 0.9853 1.0369 0.9935 0.9720 ...
             0.9999 0.994 0.969 1.0081 0.9824 0.9906 0.9776];
    f_opt = -8249.0629;
    
elseif P==177
    % 4 linear inequalities 
    % 11 nonlinear inequalities    
    Name='HS 388';
    b_L=[-70 -361 -265 -395];
    b_U=inf*ones(4,1);
    A=-[1    2  3  4   5   6   7   8   9  10 15 16 17 18 19;
        45  25 35 85  40  73  17  52  86  14 30 50 40 70 60;
        53  74 26 17  25  25  26  24  85  35 14 23 37 56 10;
        12  43 51 39  58  42  60  20  40  80 75 85 95 23 67]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    x_opt = [0.627;1.433;1.463;0.731;0.786;1.205;-1.143;1.061;-0.134;1.182;0.969;-0.845;0.481;-0.34;0.686]';
    f_opt = -5821.08;
    
elseif P==178
    % 4 linear inequalities 
    % 11 nonlinear inequalities    
    Name='HS 389';
    b_L=[-70 -361 -265 -395];
    b_U=inf*ones(4,1);
    A=-[1    2  3  4   5   6   7   8   9  10 15 16 17 18 19;
        45  25 35 85  40  73  17  52  86  14 30 50 40 70 60;
        53  74 26 17  25  25  26  24  85  35 14 23 37 56 10;
        12  43 51 39  58  42  60  20  40  80 75 85 95 23 67]; 
    c_L = [0;0;0;0;0;0;0;0;0;0;0];
    c_U = [];
    x_0 = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; %Feasible
    x_L = [];
    x_U = []; 
    x_min = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2];
    x_max = [10;10;10;10;10;10;10;10;10;10;10;10;10;10;10];
    x_opt = [0.671;1.388;1.468;0.76;0.829;1.164;-1.258;0.982;0.068;1.147;0.986;-0.888;0.565;-0.581;0.721]';
    f_opt = -5809.72;
    
elseif P==179
    % 1 nonlinear inequality
    Name='HS 394';
    b_L=[];b_U=[];A=[]; 
    c_L = 1;
    c_U = 1;
    x_0 = [2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2]; %Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
    x_opt = [0.91287;0.408268;-0.000017;-0.0000054;0.000002;-0.0000089;0.0000082;-0.000014;...
            0.000022;-0.000014;0.0000135;-0.000004;-0.000011;-0.000013;0.0000079;0.00002;...
            0.00000456;-0.000009;-0.00001;-0.0000014]';
    f_opt = 1.91667;
    
elseif P==180
    % 1 nonlinear inequality
    Name='HS 395';
    b_L=[];b_U=[];A=[]; 
    c_L = 1;
    c_U = 1;
    x_0 = [2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;...
            2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2;2]; %Not feasible
    x_L = [];
    x_U = []; 
    x_min = [-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;...
            -1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1;-1];
    x_max = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;...
            1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
    x_opt = [0.91285;0.40829;-0.0000065;-0.00000991;0.000119;-0.0000465;...
            0.0000576;-0.000048;0.0000257;0.0000117;-0.000031;0.0000087;...
            0.00002;-0.000012;-0.0000164;0.00000734;0.000017;0.0000044;...
            -0.0000059;-0.0000025;0.0000046;0.00000325;-0.00000666;-0.0000144;...
            -0.000012;-0.0000039;0.00000099;0.00000015;0.00000068;0.0000024;...
            0.0000054;0.0000027;-0.00000293;-0.000038;0.00000061;0.0000044;...
            0.0000041;0.000001455;-0.00000126;-0.000003;-0.00000386;-0.00000426;...
            -0.00000451;-0.0000045;-0.0000038;-0.00000234;-0.00000075;-0.0000000546;-0.0000011;-0.0000021]';
    f_opt = 1.91667;
end

% Define the Prob
c  = 'chs_c';
dc = 'chs_dc';

if isempty(c_L) & isempty(c_U)
    c  = [];
    dc = [];
end

% Define the Prob
Prob = conAssign('chs_f','chs_g','chs_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, A, b_L, b_U, c, dc,...
    [], ConsPattern, c_L, c_U, x_min, x_max, f_opt, x_opt);
Prob.P      = P;

if P==139 | P==156
    Prob.FUNCS.dc = '';
end
if P==152
    Prob.FUNCS.H = '';
end

% MODIFICATION LOG:
%
% 990601  mbk  First version.
% 990610       Richardo Augilar & Henrik Jönsson
% 010416  hkh  Major revision for v3.0
% 031114  ango Checked problem 167, added ConsPattern
% 040125  hkh  Corrected solution to P=79, HS 109, 2 digits switched
% 041117  med  xxx_prob removed and code added
% 060814  med  FUNCS used for callbacks instead
% 080603  med  Switched to conAssign, cleaned
% 081016  nhq  Switched all x_opt into row-vectors
% 081016  nhq  Changed f_opt for chs_prob 41, old value not correct.
% 081016  nhq  Changed f_opt for several chs_probs, old values not correct.
% 081029  nhq  Found small mistakes in several problems from chs_prob1.
% 081030  nhq  Found small mistakes in several problems from chs_prob2.
% 081031  nhq  Found small mistakes in several problems from chs_prob3.
% 081101  nhq  Found small mistakes in several problems from chs_prob4.
