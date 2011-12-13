% uc_prob: Defines unconstrained optimization problems (with simple bounds)
%
% function [probList, Prob] = uc_prob(P);
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
% Written June 1, 1999.   Last modified Aug 4, 2008.

function [probList, Prob] = uc_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Rosenbrocks banana'...
    ,'Exponential'...
    ,'Polynomial'...
    ,'Spiral'...
    ,'Fletcher Q.2.3'...
    ,'Fletcher Q.3.9 Quadratic'...
    ,'AG-1. Symmetric exp.func. 2 min'...
    ,'AG-2. Deep hole in origo'...
    ,'RE-1. Marsden/Tromba fig. 2.1.17'...
    ,'RE-2. Netlib/trig, (n=2)'...
    ,'RE-3. Structured exponential, (p=4)'...
    ,'Fletcher Q.2.4'...
    ,'Fletcher Q.2.5'...
    ,'Fletcher Q.2.2'...
    ,'Fletcher Q.2.7 Quadratic'...
    ,'Fletcher Q.2.6'...
    ,'Fletcher Q.3.3'...
    ,'3rd degree'...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

HessPattern = []; pSepFunc = []; uP = [];
f_opt = []; x_opt = []; f_Low = [];

if P == 1
    Name='Rosenbrocks banana';
    x_0=[-1.2 1]';
    x_L=[-10;-10];
    x_U=[2;2];
    x_opt=[1 1];
    f_opt=0;
    f_Low=0;
    x_max=[ 1.3  1.3];
    x_min=[-1.1  -0.2];
    pSepFunc=2;
elseif P == 2
    Name='Exponential Gill-Murray-Wright';
    x_0=[1 -0.5]';
    x_L=[-100;-100];
    x_U=[100;100];
    % Optimal point at (0.5, -1). Saddle at (-1.5,1), Min at (-inf,*)
    x_opt=[0.5 -1 0;-1.5 1 1];
    f_opt=[0;1.785];
    f_Low=-10;
    x_max=[ 0.6  1.1];
    x_min=[-1.6 -1.5];
elseif P == 3
    Name='Polynomial Function';
    x_0=[1.5 1 0.5 ]';
    x_L=[1E-14;1E-14;1E-14]';
    x_U=[1000;1000;1000]';
    x_opt=[1.5375 1.1138 0.55692];
    f_Low=0;
    x_max=[ 2  2  2];
    x_min=[0.001 0.001 0.001];
elseif P == 4
    Name='Spiral function';
    uP(2,1) = 0.4;
    uP(1,1) = 1;
    x_0=[-0.5 1]';
    x_L=[-100;-100];
    x_U=[100;100];
    f_Low=-100;
    x_max=[8 8];
    x_min=[-8 -8];
elseif P == 5
    Name='Fletcher Q.2.3';
    x_0=[0.5 0.5]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[ 1 0 0;0 0 1;0 -1 1;-1 -1 2]; % 4 Stationary points, one per row
    f_Low=-100;
    x_max=[ 1.5  0.5];
    x_min=[-1.5 -1.5];
elseif P == 6
    % Fletcher Q.3.9. QP. BFGS converges in 2 steps with exact line search
    Name='Fletcher Q.3.9 Quadratic';
    x_0=[0 0]';
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=-100;
    x_max=[ 3  3];
    x_min=[-1 -1];
elseif P == 7
    Name='AG-1. Symmetric exponential';
    x_0=[1 -1]';
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=0;
    x_max=[ 1  1];
    x_min=[-1 -1];
elseif P == 8
    Name='AG-2. Deep hole in origo';
    x_0=[1 1]';
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=-10;
    x_max=[ 1  1];
    x_min=[-1 -1];
elseif P == 9
    %Marsden - Tromba fig. 2.1.17
    Name='RE-1 Marsden/Tromba, max';
    uP = -1; %or uP = 1
    f_opt = -3; %or f_opt = 0
    x_0=[0.5 0.5]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[];
    x_max=[ 2  2];
    x_min=[-2 -2];
elseif P == 10
    %netlib/trig for n=2
    Name='RE-2 Netlib/trig';
    x_0=[1 1]';
    x_L=[-10;-10];
    x_U=[10;10];
    f_Low=-1;
    x_max=[ 6.5  6.5];
    x_min=[-6.5 -6.5];
elseif P == 11
    Name='Structured Exponential Gill et.al';
    x_0=[1 -0.5]';
    x_L=[-100;-100];
    x_U=[100;100];
    x_opt=[0.5 -1]; % Optimal point. Saddle at (-1.5,1), Min at (-inf,*)
    f_opt=0;
    f_Low=-10;
    x_max=[ 0.6  1.1];
    x_min=[-1.6 -1.5];
    uP(1)=0;
elseif P == 12
    Name='Fletcher Q.2.4';
    x_0=[1 0]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[0 0 1];	% (0,0) is a saddle, no local min/max
    f_opt=0;
    f_Low=-100;
    x_max=[ 2  2];
    x_min=[-2 -2];
elseif P == 13
    Name='Fletcher Q.2.5';
    x_0=[-1 -1.5]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[0 0];
    f_opt=0;
    f_Low=-10;
    x_max=[ 2  2];
    x_min=[-2 -2];
elseif P == 14
    Name='Fletcher Q.2.2';
    x_0=[1 1]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[0.69588;-1.3479];
    f_opt=-.58245;
    f_Low=-100;
    x_max=[2 2];
    x_min=[-2 -2];
elseif P == 15
    Name='Fletcher Q.2.7 Quadratic';
    x_0=[0 0]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[];
    f_opt=[];
    f_Low=-100;
    x_max=[3 3];
    x_min=[-3 -3];
elseif P == 16
    Name='Fletcher Q.2.6';
    x_0=[0 0]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[-1 1 0;-1 0 0;-1 2 0];  % Minimum along a line x(1)==1
    f_opt=-4;
    f_Low=-100;
    x_max=[1 1];
    x_min=[-1 -1];
elseif P == 17
    Name='Fletcher Q.3.3';
    x_0=[1 1]';
    x_L=[-10;-10];
    x_U=[10;10];
    x_opt=[0 0 0; 0 1 1];
    f_opt=0;
    f_Low=-100;
    x_max=[2 2];
    x_min=[-2 -2];
elseif P == 18
    Name='Nash-Sofer page 307 #3';
    x_0=[1 1 ]';
    x_L=[-10;-10]';
    x_U=[100;100]';
    x_opt=[];
    f_opt=[];
    f_Low=-1000;
    x_max=[ 2 2];
    x_min=[-2 2];
elseif P == 19
    Name='3rd degree';
    % (0,0) minimum, (1,2) saddle, (5/3,0) maximum
    x_0=[1 1 ]';
    x_L=[-10;-10]';
    x_U=[100;100]';
    x_opt=[0 0];
    f_opt=[0];
    f_Low=-1000;
    x_max=[ 2 2];
    x_min=[-2 2];
else
    error('uc_prob: Illegal problem number');
end

% Define the Prob
Prob = conAssign('uc_f','uc_g','uc_H', HessPattern, x_L,...
    x_U, Name, x_0, pSepFunc, f_Low, [], [], [], [], [],...
    [], [], [], [], x_min, x_max, f_opt, x_opt);
Prob.P  = P;
Prob.uP = uP;

% MODIFICATION LOG:
%
% 980826  hkh  Defining probType before call to ProbVarDef.
% 980820  hkh  Change name f_min to f_Low. Delete dummy problems HKH-2-HKH-5
% 981005  mbk  b_L=[]; b_U=[]; c_L=[]; c_U=[]; for all problems.
% 981006  hkh  Added call to checkuP
% 981010  hkh  Changed to use tomFiles for name definitions
% 981022  hkh  Set Name=[] if P is illegal
% 981027  hkh  Check which P in Prob.P, not just on if nonempty
% 981102  hkh  Shorter names for two functions
% 990909  hkh  Add HessPattern
% 011030  hkh  Change problem names, to avoid too long names
% 041117  med  xxx_prob removed and code added
% 080603  med  Switched to conAssign, cleaned
% 080804  hkh  Added simple 3rd degree problem as #19
