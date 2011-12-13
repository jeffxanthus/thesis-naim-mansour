% gp_prob: Defines Geometric Programming problems
%
% function [probList, Prob] = gp_prob(P);
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

function [probList, Prob] = gp_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Dembo78' ...
    ,'HANDA' ...
    ,'P01' ...
    ,'P02' ...
    ,'P03' ...
    ,'P04' ...
    ,'P05' ...
    ,'P06' ...
    ,'P07' ...
    ,'P08' ...
    ,'P09' ...
    ,'P10A' ...
    ,'P11' ...
    ,'P12'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return;
end

if P==1
    % Dembo78. COPL_GP test problem
    Name='Dembo78';
    GP.nterm = [2;2];
    GP.coef = [1;1;1/4;1];
    GP.A = sparse([1 -1 0.5 0; 1 -1 0 1]');
    x_L = zeros(2,1);
    x_U = x_L + inf;
    x_0=0.25*ones(2,1);	% Starting values for the optimization
    x_min=zeros(2,1);
    x_max=2*ones(2,1);
    x_opt=[4.219748195 .236980965];
    f_opt=2;
elseif P==2
    % handa.dat
    Name='HANDA';
    GP.nterm = [6 2 3 2]';
    GP.coef = ...
        [.1e1 .25e1 .3e1 .7e1 .9 .1 .4 .1e1 .2e1 .1e1 .1e1 .5 .3e1]';
    GP.A = sparse([ 2   0  2   0  0  0  0  0  0  0;...
        2   0  0   1 -1  0  0  0  0  0;...
        1  -1  0   1  1  0  0  0  0  0;...
        -6  -1  0  -3  0  0  0  0  0  0;...
        1/3   0  0   0  0  1  0  1  0  0;...
        0  -1  0   0  0  0  0  0  1  0;...
        -1   1  0   0  0  0  0  0  0  0;...
        1   0 -1   0  0  0  0  0  0  0;...
        -3   0 -1/3 0  0  1 -3  0  0  0;...
        -1   0  0   0  0 -1  1 -1  0  0;...
        2  -1  0   0  0  0  0  3  0  0;...
        -1   0 -1   0  0  0  0  0  0  1;...
        0  -1 -1   0  0  0  0  0 -1 -1.5]);
    x_L = zeros(10,1);
    x_U = x_L + inf;
    x_0=0.25*ones(10,1);	% Starting values for the optimization
    x_min=zeros(10,1);
    x_max=15*ones(10,1);
    x_opt=[1.1656891;1.1690029;1.9465039;1.0295003;1.0649260;11.224165;...
        3.6643909;.40449645;.73360368;2.7226570];
    f_opt=.18262629e2;
elseif P==3
    % P1
    Name='P01';
    GP.nterm = [6;3];
    GP.coef = [.5e1;.5e5;.2e2;.72e5;.1e2;.144e6;.4e1;.32e2;.12e3];
    GP.A = sparse([ 1  -1  0   0  0  0 -1  0  0;...
        0   0  1  -1  0  0  0 -1  0;...
        0   0  0   0  1 -1  0  0 -1])';
    x_L = zeros(3,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(3,1);
    x_max=250*ones(3,1);
    x_opt=[108.66933;85.147491;204.29971];
    f_opt=.62998434e4;
elseif P==4
    % P2
    Name='P02';
    GP.nterm = [2; 1; 2];
    GP.coef = [1;25;1;5;1];
    GP.A = sparse([ -1   1  1   1  0;...
        -1  -1  1   0  1;...
        1   0  0  -1 -1])';
    x_L = zeros(3,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(3,1);
    x_max=2*ones(3,1);
    x_opt=[.2418710273;4.13443483;5.343789967];
    f_opt=6.806329814;
elseif P==5
    % P3
    Name='P03';
    GP.nterm = [3; 1; 2];
    GP.coef = [.49;200;.73;.274e7;.6667e-6;1];
    GP.A = sparse([ 1   1  1  -1  0 0;...
        .6667 1 .6667 -1  1 0;...
        3 0 3 -1  1 0;...
        0 0 1  -1  1 1])';
    x_L = zeros(4,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(4,1);
    x_max=2*ones(4,1);
    x_opt=[15.153905;16513.899;12.449836;.87945302];
    f_opt=.71523492e8;
elseif P==6
    % P4
    Name='P04';
    GP.nterm = [1; 3; 3; 3];
    GP.coef = [1;2;1;3;1;3;2;1;1;1];
    GP.A = sparse([ -1  1  0  0  1  0  0  1  0  0;...
        -1  0  1  0  0  1  0  0  1  0;...
        -1  0  0  1  0  0  1  0  0  1])';
    x_L = zeros(3,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(3,1);
    x_max=2*ones(3,1);
    x_opt=[.21635386;.17376484;.13117580];
    f_opt=.20277748e3;
elseif P==7
    % P5
    Name='P05';
    GP.nterm = [5; 2; 1; 2; 1; 1; 3; 3; 1];
    W = 240;
    GP.coef = [81.833*W;216.25*W;216.25*W;8.167*W;1.402*W;.87;10.11;1.58;8.728;71.273;...
        0.008;0.0018;2.7065e-7*W^2;4.6083e-4*W^2;5.8371e-3*W^2;.9877e-8*W^3;...
        1.5687e-5*W^3;1.794e-4*W^3;.112];
    GP.A = sparse([ 1  0  0  0  0 -1 -2 -1 -2 -3  1  1  0  1  0  0  1  0  0;...
        0  1  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0;...
        0  0  1  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0;...
        0  0  0  1  0  0  0  0  0  0  0  0 -1 -2 -2 -2 -3 -3  1;...
        0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 -1])';
    x_L = zeros(5,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(5,1);
    x_max=2*ones(5,1);
    x_opt=[4.8421699;.38737356e-1;.87159052e-2;21.612075;2.4205524];
    f_opt=.14073723e6;
elseif P==8
    % P6
    Name='P06';
    GP.nterm = [1;10];
    GP.coef = [10e16;.4411;28.46;616;.03703;710.7;.3225;2.93;.04471;.3796;4.289];
    GP.A = sparse([ -2  1  2  2  0  0  1  0  0  0  1;...
        -1  0  0  1  0  0  0  1  1  2  1;...
        -1  0  0  0  1  2  1  1  0  0  0])';
    x_L = zeros(3,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(3,1);
    x_max=2*ones(3,1);
    x_opt=[.56070739e-1;.24731523;.20447030e-1];
    f_opt=.62899378e21;
elseif P==9
    % P7
    Name='P07';
    GP.nterm = [4; 3; 3; 4; 4];
    GP.coef = [10;15;20;25;.5;.7;.2;1.3;.8;3.1;2;.1;1;.65;.2;.3;.4;.5];
    GP.A = sparse([ 1 -1 -2  2 .5  3  0 -.5  0  -1    1   0  -1   0 -2 .5 -3  0;...
        -1 -2  1  2  0  1 -1   1  0  .5    0   1   1  -2  1  2 -2  0;...
        0  1  0 -1 -1 -2  1  -1  1   0 -1.5 -0.5 0.5  1  0  1  1 -2;...
        2  1 -1  0  0  0 -0.5 0 -1  -2    0   0   0   0 -1 .3333 0 1;...
        0 -1 -2 .5  0  0  0  -1 -1  -1    1   1   1   1 .5 -.6667 1 0;...
        -3  0  1 -2 -2  1 .6667 1 2 .3333 -1 -1 0 -1 0 0 0 0;...
        -.25 -.5 0 1 1 .5 .25 0 0 0 .3333 -.5 0 1 .3333 .25 .75 .5])';
    x_L = zeros(7,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(7,1);
    x_max=2*ones(7,1);
    x_opt=[2.8562502;.61086850;2.1510043;4.7123256;.99933060;1.3472654;.31652931e-1];
    f_opt=.18099866e4;
elseif P==10
    % P8
    Name='P08';
    GP.nterm = [4; 3; 3; 4; 4];
    GP.coef = [10;15;20;25;.5;.7;.2;1.3;.8;3.1;2;.1;1;.65;.2;.3;.4;.5];
    GP.A = sparse([ 1 -1 -2  2 .5  3  0 -.5  0  -1    1   0  -1   0 -2 .5 -3  0;...
        -1 -2  1  2  0  1 -1   1  0  .5    0   1   1  -2  1  2 -2  0;...
        0  1  0 -1 -1 -2  1  -1  1   0 -1.5 -0.5 0.5  1  0  1  1 -2;...
        2  1 -1  0  0  0 -0.5 0 -1  -2    0   0   0   0 -1 .3333 0 1;...
        0 -1 -2 .5  0  0  0  -1 -1  -1    1   1   1   1 .5 -.6667 1 0;...
        -3  0  1 -2 -2  1 .6667 1 2 .3333 -1 -1 0 -1 0 0 0 0;...
        .125 -.5 0 1 1 .5 .25 0 0 0 .3333 -.5 0 1 .3333 .25 .75 .5])';
    x_L = zeros(7,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(7,1);
    x_max=2*ones(7,1);
    x_opt=[3.8969138;.80946213;2.6649111;4.3000107;.85335350;1.0949624;.27313093e-1];
    f_opt=.91199012e3;
elseif P==11
    % P9
    Name='P09';
    GP.nterm = [4; 3; 3; 4; 4];
    GP.coef = [10;15;20;25;.5;.7;.2;1.3;.8;3.1;2;.1;1;.65;.2;.3;.4;.5];
    GP.A = sparse([ 1 -1 -2  2 .5  3  0 -.5  0  -1    1   0  -1   0 -2 .5 -3  0;...
        -1 -2  1  2  0  1 -1   1  0  .5    0   1   1  -2  1  2 -2  0;...
        0  1  0 -1 -1 -2  1  -1  1   0 -1.5 -0.5 0.5  1  0  1  1 -2;...
        2  1 -1  0  0  0 -0.5 0 -1  -2    0   0   0   0 -1 .3333 0 1;...
        0 -1 -2 .5  0  0  0  -1 -1  -1    1   1   1   1 .5 -.6667 1 0;...
        -3  0  1 -2 -2  1 .6667 1 2 .3333 -1 -1 0 -1 0 0 0 0;...
        .5 -.5 0 1 1 .5 .25 0 0 0 .3333 -.5 0 1 .3333 .25 .75 .5])';
    x_L = zeros(7,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(7,1);
    x_max=2*ones(7,1);
    x_opt=[4.3938838;.85445219;2.8431382;3.4002319;.72287845;.87041068;.24643648e-1];
    f_opt=.54373613e3;
elseif P==12
    % P10A
    Name='P10A';
    GP.nterm = [3; 3; 1; 1; 1; 1; 1; 1];
    GP.coef = [2;5;4.7;7.2;.5;.2;10;.6;6.2;3.1;3.7;.3];
    GP.A = sparse([ .9  0  0 -3.8 0   0 2.3  0    0  1.6   0    0;...
        -1.5 0  0  2.2 0   0 1.7  0    0  .4    0    0;...
        -3  0  0  4.3 0   0 4.5 -2.1  0 -3.8   0    0;...
        0 -.3 0   0 -.7  0  0  -2.1  0   0   5.4   0;...
        0 2.6 0   0 -1.6 0  0  .4    0   0   1.3   0;...
        0  0 -1.8 0  0  4.3 0   0   4.5  0    0   -1.1;...
        0  0 -.5  0  0 -1.9 0   0  -2.7  0    0    7.3;...
        0  0   1  0  0  8.5 0   0   -.6  0    0   -5.6])';
    x_L = zeros(8,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(8,1);
    x_max=2*ones(8,1);
    x_opt=[.96889019;.19895172;1.1212709;.78441002;1.0022438;.70103910;1.0941521;.97245608];
    f_opt=.29229484e2;
elseif P==13
    % P11
    Name='P11';
    GP.nterm = [8;3;7;7;7;6;5;5];
    GP.coef = [10;20;30;100;5;50;25;10;.1;.05;.15;.1;.05;.05;.15;.1;.1;.2;...
        .1;.1;.1;.1;.1;.1;.1;.02;.02;.02;.02;.02;.02;.02;.01;.01;...
        .01;.01;.01;.01;.1;.1;.1;.1;.1;.02;.02;.02;.02;.02];
    GP.A = sparse([ 1 -1 0 -1   2    0  0  0 2  0  0  1  1  0 -.5 0 0 1 1 0 0 0 0 0 0 2 1 1 1 1 1 1  1  0  0  0  0  0  1  0  0  0  0 -.5   0   0   0   0;...
        -1  0 1 -1   2    0  0  0 2  0  0  0 -1  2 -.3 0 0 1 0 1 0 0 0 0 0 0 1 0 0 0 0 0 -1  1  0  0  0  0  0  1  0  0  0   0 -.5   0   0   0;...
        -1  0 1 -1   1  -.5  2 .5 1  0  0  0 -1  2   1 0 0 1 0 0 1 0 0 0 0 0 0 1 0 0 0 0  0 -1  1  0  0  0 -2  0  1  0  0   1   0 -.5   0   0;...
        1 -1 1 -1   0  -.5  2 .5 0  1  0  1  0 -1   0 0 2 0 0 0 0 1 0 0 0 0 0 0 1 0 0 0  0  0 -1  1  0  0  0 -2  0  1  0   0   1   0 -.5   0;...
        0  1 0 -1   1  -.5 -1  1 0 .5  0  0  1  0  .5 1 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0 0  0  0  0 -1  1  0  0  0 -2  0  1   0   0   1   0 -.6;...
        0  1 0 -1 1.5    0 -1  1 0  0 .5  0  1  0   0 1 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1 0  0  0  0  0 -1  1  0  0  0 -2  0   0   0   0   1   0;...
        0  0 0 -1   2    0 -1  1 0  0 .5  1 .5  0   0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 1  0  0  0  0  0 -1  0  0  0  0 -2   0   0   0   0  -1])';
    x_L = zeros(7,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(7,1);
    x_max=2*ones(7,1);
    x_opt=[1.3417762;.99321692;.87041920;.92347603;3.1489302;.40380572;1.5485695];
    f_opt=.17847793e3;
elseif P==14
    % P12
    Name='P12';
    GP.nterm = [1;4;12;14];
    GP.coef = [100000;.05367;.02186;.09773;.006694;1e-6;1e-5;1e-6;1e-10;1e-8;1e-2;1e-4;...
        .109;1.611e-4;1e-23;1.93e-6;1e-3;1e-6;1e-5;1e-6;1e-9;1e-9;1e-3;1e-3;.109;...
        1.611e-5;1e-23;1.93e-8;1e-5;1.118e-4;1e-4];
    GP.A = sparse([ -.0013 1 0 0 0 1 0 0 0  0  0 0 0  0 0  0  0 1 0 0 0  0 0 0 0 0 0  0 0 1 0;...
        -.0023 0 1 0 0 0 1 0 0  0  0 0 0  1 1  1  0 0 1 0 0  0 0 0 0 1 1  1 0 0 0;...
        -.0025 0 0 1 0 0 0 1 0  0  0 0 0  0 0  0  0 0 0 1 0  0 0 0 0 0 0  0 0 0 0;...
        -4.67  0 0 0 1 0 0 0 1  0  0 0 1  0 1 -1  0 0 0 0 1  0 0 0 1 0 1 -1 0 0 0;...
        -4.672 0 0 0 1 0 0 0 0  1  0 0 1  1 1  1  0 0 0 0 0  1 0 0 1 1 1  1 0 0 0;...
        -.0081 0 0 0 0 0 0 0 0  0  1 0 0  0 0  0  0 0 0 0 0  0 1 0 0 0 0  0 0 0 0;...
        -.0081 0 0 0 0 0 0 0 0  0  0 1 0  0 0  0  0 0 0 0 0  0 0 0 0 0 0  0 0 0 0;...
        -.005  0 0 0 0 0 0 0 0  0  0 0 0  0 0  0  0 0 0 0 0  0 0 1 0 0 0  0 0 0 0;...
        -.0009 0 0 0 0 0 0 0 0  0  0 0 0  0 0  0  0 0 0 0 0  0 0 0 0 0 0  0 1 1 0;...
        -.0009 0 0 0 0 0 0 0 0  0  0 0 0  0 0  0  1 0 0 0 0  0 0 0 0 0 0  0 0 0 0;...
        -.0012 0 0 0 0 0 0 0 0  0  0 0 0  0 0  0  0 0 0 0 0  0 0 0 0 0 0  0 0 0 1;...
        0 0 0 0 0 0 0 0 0 -1 -1 1 0 -1 0 -2 -1 0 0 0 0 -1 0 0 0 0 0  0 0 0 0])';
    x_L = zeros(12,1);
    x_U = x_L + inf;
    x_0=0.25*ones(length(x_L),1);	% Starting values for the optimization
    x_min=zeros(12,1);
    x_max=2*ones(12,1);
    x_opt=[2.4025977;2.5646752;7.7143203;1.1809676;7.7248836;1.3009379;4.2562875;2.7820972;1.8025610;2.0493149;6.7021688;6.5437335];
    f_opt=.31697051e1;
else
    error('gp_prob: Illegal problem number');
end

Prob = gpAssign(GP.nterm, GP.coef, GP.A, Name);
Prob.x_0 = x_0;
Prob.x_L = x_L;
Prob.x_U = x_U;
Prob.x_min = x_min;
Prob.x_max = x_max;
Prob.x_opt = x_opt;
Prob.f_opt = f_opt;

% MODIFICATION LOG:
%
% 050509  med  Created
% 050510  med  Another problem added
% 050519  med  More problems added (P01 to P12)
% 050520  med  Format changed, first column removed, GP.A now transposed
% 050608  med  Added gp_c for nonlinear constraints (primal)
% 080603  med  Switched to gpAssign, cleaned