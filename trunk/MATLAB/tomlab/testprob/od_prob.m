% od_prob: Defines parameter estimation in ODEs as constrained
%          nonlinear least squares problems
%
% function [probList, Prob] = od_prob(P);
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

function [probList, Prob] = od_prob(P, varargin)

if nargin < 1
    P=[];
end

probList=str2mat(...
    'Arporn Reactor' ...
    ); % MAKE COPIES OF THE PREVIOUS ROW AND CHANGE TO NEW NAMES

if isempty(P)
    return
end

f_Low = []; xName = []; x_opt = []; f_opt = [];
cName = []; A = []; b_L = []; b_U = []; c_L = []; c_U = [];
y = []; SepAlg = []; weightType = []; weightY = [];
JacPattern = []; ConsPattern = []; pSepFunc = []; uP = []; uPName = [];

param = []; Z = [];

if P==1
    Name='Arporn Reactor';
    % Arporn Amornchai reactor fitting problem, 3 short data series
    x_0 = [0.8 100 0.2 100]';
    x_L = [0 0 0 0]';
    x_U = [50 30000 100 30000]';
    x_min = x_L;
    x_max = x_L;
    E = [
        362.96 365.03 366.11
        366.14 369.03 369.78
        365.57 369.04 369.54
        369.15 372.66 373.46
        390.49 392.49 392.72
        385.82 389.76 390.54 ];
    t = [0.5 1.9 3.3 4.7 6.1 7.5]';
    Eeq = [6 6 6];
    EWht      = [1 1 1];
    resWeight = [];
    odePattern = [
        0 0 0 0
        1 1 1 1
        1 1 0 0
        1 1 1 1
        0 0 1 1
        1 1 1 1];
    Y0 = [1 0.08 1.15 1.35 2.76 273+90];
    odeT_s  =   0;
    odeT_e  = 7.9;
    odeH_s  = [];
    void        = 0.4;
    Lvoid       = 0.33;
    Gvoid       = 0.07;
    ug          = 178.77/Gvoid;   % m/hr
    ul          = 34.87/Lvoid;    % m/hr
    kial        = 1.97e-3*3600;   % 1/hr         %  8.79e-4*3600
    Cpg         = 29.3e3/2;       % J/(kg K)
    Cpl         = 2.17e3/84.4;    % J/(kg K)
    rhog        = 2;              % kg/m3
    rhol        = 820;            % kg/m3
    param.H     = -125.52e6;      % J/kmol
    param.R     = 8.314;
    param.valG  = -(1-void)*kial/(Gvoid*ug);
    param.valL1 = -(1-void)/(Lvoid*ul);
    param.valL2 = (1-void)*kial/(Lvoid*ul);
    param.valL3 = (1-void)/(Lvoid*ul);
    param.valT  = (1-void)/(ug*Gvoid*rhog*Cpg + ul*Lvoid*rhol*Cpl);
else
    error('od_prob: Illegal problem number');
end

Prob=clsAssign();

Prob = odeFitAssign(Prob, '', Y0, E, Eeq, x_L, x_U, Name, ...
    t, x_0, odePattern, EWht, resWeight,...
    odeT_s, odeT_e, odeH_s);

% Send constant parameters in Prob.ODE.user
Prob.ODE.param = param;
Prob.ODE.Z     = Z;

% Defined ODE function routine
Prob.ODE.f = 'odFunc_f';
Prob.ODE.J = '';

% MODIFICATION LOG:
%
% 050418  hkh  Written, based on cls_prob
% 050503  hkh  Change to probType ODE
% 080603  med  Switched to conAssign, cleaned