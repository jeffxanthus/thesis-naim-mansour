% goals_prob: Defines unconstrained and constrained Multi criterium problems
%
% function [probList, Prob] = goals_prob(P);
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

function [probList, Prob] = goals_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'EASY-TP269'...
    ,'EASY-TP13'...
    ,'EASY-TP46'...
    ,'EASY-TP48'...
    ,'EASY-TP354'...
    ,'EASY-TP355'...
    ,'EASY-TP372'...
    ,'EASY-TP373'...
    ,'EASY-OPT_KINX'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

JacPattern = []; t = []; weightType = []; weightY = [];
c_L = []; c_U = []; SepAlg = []; f_Low = [];
ConsPattern = []; A = []; b_L = []; b_U = [];

if P==1
    Name='EASY-TP269';
    %Constrained least squares problem with four linear terms
    %Hock W., Schittkowski K. Test examples for Nonlinear Programming
    %m=3; n=5;
    x_0=[2 2 2 2 2]';
    x_opt = [-0.76744186 0.25581395 0.62790698 -0.11627907 0.25581395]'; % estimated from eASY-FIT
    f_opt =0.20465116e+01; % Estimated from EASY-FIT
    x_L = [-1e3;-1e3;-1e3;-1e3;-1e3];
    x_U = [1e3;1e3;1e3;1e3;1e5];
    x_min =[];
    x_max =[];
    A= [1 3 0 0 0;0 0 1 1 -2;0 1 0 0 -1];
    b_L = zeros(3,1);
    b_U = zeros(3,1);
    y = zeros(4,1);
elseif P==2
    Name='EASY-TP13';
    %Hock W., Schittkowski K. (1981). Test examples for Nonlinear Programming
    x_0=[0 0]';
    x_opt = [0.99229265 0]'; % Estimated from EASY-FIT
    f_opt =0.50773705; % Estimated from EASY-FIT
    x_L = [0;0];
    x_U = [1.5;1.5];
    x_min =[];
    x_max =[];
    c_L = 0;
    c_U = 0;
    y = zeros(2,1);
elseif P==3
    Name='EASY-TP46';
    % Equality constrained academic test problem
    % Hock W., Schittkowski K. (1981):
    % Test examples for Nonlinear Programming
    x_0=[1.007366   1.007399   1.000027   1.015641   1.046740]';
    x_opt = [0.99930297 0.99930297 1.0000000 1.0003484 0.99930284]';
    f_opt =0.43388180e-14;% Estimated from EASY-FIT
    x_L = zeros(5,1);
    x_U = 1e5*ones(5,1);
    x_min =[];
    x_max =[];
    c_L = [1 2]';
    c_U = [1 2]';
    y = zeros(4,1);
elseif P==4
    Name='EASY-TP48';
    %Equality constrained academic test problem
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    x_0 = [1.000000   1.000000   1.000000   1.000000   1.000000]';
    x_opt = [0.99998200 -3.0000167 -3.0000223 5.0000229 5.0000341]'; % Estimated from EASY-FIT
    f_opt =0.47988882e-9;
    x_L = -1e3*(ones(5,1));
    x_U = 1e5*(ones(5,1));
    x_min =[];
    x_max =[];
    A = [1 1 1 1 1;0 0 1 -2 -2];
    b_L = [5;-3];
    b_U = [5;-3];
    y = zeros(3,1);
elseif P==5
    Name='EASY-TP354';
    %Constrained least squares problem, four quadratic terms
    %Hock W., Schittkowski K.Test examples for Nonlinear Programming
    x_0 = [3 -1 0 1]';
    x_opt = [0.50331793 -0.45561167e-1 0.23581898 0.30642425]'; % Estimated from EASY-FIT
    f_opt = 0.56891924e-01;
    x_L = -1e3*ones(4,1);
    x_U = 1e5*ones(4,1);
    x_min =[];
    x_max =[];
    A= [1 1 1 1] ;
    b_L =1;
    b_U =1;
    y = zeros(4,1);
elseif P==6
    Name='EASY-TP355';
    %Constrained least squares problem, four quadratic terms and local solutions
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
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
    y = zeros(2,1);
elseif P==7
    Name='EASY-TP372';
    %Least squares problem, twelve inequality constraints
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    x_0=[3.0000e2 -1.0000e2 -2e-1 -1.2700e2 -1.5100e2 3.7900e2 4.2100e2 4.600e2 4.2600e2]';
    x_opt = [523.30554 -156.95784 -0.19966457 29.608010 86.615555 47.326700 26.235597 22.915985 39.470737]'; % Estimated from EASY-FIT
    f_opt = 0.18620249e-30;  % Estimated from EASY-FIT Estimated from MINOS/snopt fopt_= 6695.046559741104800000
    x_L =-1e4*ones(9,1);
    x_U = 1e5*ones(9,1);
    x_min =[];
    x_max =[];
    c_L =[127;151;379;421;460;426;-127;-151;-379;-421;-460;-426];
    c_U =[];
    y = zeros(6,1);
elseif P==8
    Name='EASY-TP373';
    %Least squares problem, six inequality constraints
    %Hock W., Schittkowski K. (1981):
    %Test examples for Nonlinear Programming
    x_0=[3.0000e2 -1.0000e2 -2e-1 -1.2700e2 -1.5100e2 3.7900e2 4.2100e2 4.600e2 4.2600e2]';
    x_opt = [523.30554 -156.95784 -0.19966457 29.608010 -86.615555 47.326700 26.235597 22.915985 -39.470737]'; % Estimated from EASY-FIT
    f_opt = 0.18620249e-30;
    x_L =-1e4*ones(9,1);
    x_U = 1e5*ones(9,1);
    x_min =[];
    x_max =[];
    c_L = [127;151;379;421;460;426];
    c_U = [];
    y = zeros(6,1);
elseif P==9
    Name='EASY-OPT_KINX';
    %Schittkowski
    %Linear kinetics with variable switching times (optimal control problem)
    t = zeros(1,100);
    for i=1:100
        t(i) = i;
    end
    y = 4*ones(100,1);
    x_0 = [5.00e1 5.00e1 5.00e0 5.00e0 1.00e1 3.00e1]';
    x_opt = [39.736101 41.531040 17.706125 2.0840358 10.627758 34.767734]';
    f_opt = 0.79185898e1;
    x_L = [1 1 0 0 5 5]';
    x_U = [1e5 1e5 5e5 1e4 5e1 1e5]';
    x_min =[];
    x_max =[];
    A = [0 0 0 -1 1 0;0 0 0 0 -1 1];
    b_L = [];
    b_U = [0;0];
else
    error('goals_prob: Illegal problem number');
end

% Define the Prob
c  = 'goals_c';
dc = 'goals_dc';
if isempty(c_L) & isempty(c_U)
    c  = [];
    dc  = [];
end

Prob = clsAssign('goals_r', 'goals_J', JacPattern, x_L, x_U, Name, x_0, ...
    y, t, weightType, weightY, SepAlg, f_Low, ...
    A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
    x_min, x_max, f_opt, x_opt);
Prob.P = P;

if P == 9
    Prob.FUNCS.J = [];
end

% MODIFICATION LOG:
%
% 040517  med  Created.
% 041117  med  xxx_prob removed and code added
% 060814  med  FUNCS used for callbacks instead
% 080603  med  Switched to clsAssign, cleaned