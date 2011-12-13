% odeAssign is a direct way of setting up a Ordinary Differential Equation (ODE)
% in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% function Prob = odeAssign(f, tInit, Y0, tWant, tStop, InitStep, J, Name, Prob)
%
% It is then possible to solve the ODE using any ODE solver
%
% INPUT (At least three input parameters needed)
%   f        Name of m-file that defines y'=f(y,t)
%   tInit    Initial value of the independent variable
%   Y0       Vector of initial solutions
%            nEq (length(Y0): Number of equations in the ODE system.
% .....................................................................
%   tWant    Vector of points where a solution is requiered
%
%   tStop    If tWant is given, tStop is a limit that the solver is not
%            allowed to integrate past.
%            If tWant is not given, the solver integrates from tInit to tStop
%
%   InitStep Initial step length
%
%   J        Name of m-file that defines the Jacobian of the system
%
%   Name     Name of the ODE problem (equation)
%
%   Prob     A TOMLAB Prob structure, ODE information is set in Prob.ODE
%            If empty, an [] Prob structure is defined by a call to ProbDef(1)

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.8.0$
% Written Apr 21, 2005.  Last modified Jul 6, 2005.

function Prob = odeAssign(f, tInit, Y0, tWant, tStop, InitStep, J, Name, Prob)

if nargin < 9 
   Prob = [];
   if nargin < 8 
     Name = '';
     if nargin < 7
       J = [];
       if nargin < 6
          InitStep = [];
          if nargin < 5
             tStop = [];
             if nargin < 4
                tWant = [];
                if nargin < 3
                   error(['odeAssign requires at least three parameters, '...
                          ,'f, tInit, Y0']);
end, end, end, end, end, end, end

if isempty(Prob) | ~isstruct(Prob)
   Prob     = ProbDef(1);
   Prob.ODE = struct('P',double(1), 'absTol',[], 'relTol',[], 'f',[], 'J',[],...
                     'InitStep',[],'Name',[],'tInit',[],'tStop',[],'tWant',[],...
                     'Y0',[],'LSODE',[],'RKSUITE',[],'ML',[],'CHECK',1);
elseif ~isfield(Prob,'ODE')    
   Prob.ODE = struct('P',double(1), 'absTol',[], 'relTol',[], 'f',[], 'J',[],...
                     'InitStep',[],'Name',[],'tInit',[],'tStop',[],'tWant',[],...
                     'Y0',[],'LSODE',[],'RKSUITE',[],'ML',[],'CHECK',1);
end

Prob.ODE.f        = deblank(f);
Prob.ODE.J        = deblank(J);
Prob.ODE.InitStep = InitStep;
Prob.ODE.Name     = deblank(Name);
Prob.ODE.tInit    = tInit;
Prob.ODE.tStop    = tStop;
Prob.ODE.tWant    = tWant;
Prob.ODE.Y0       = Y0;

% MODIFICATION LOG
%
% 050407  med  Help updated
% 050412  bjo  Moved f,J from Prob.FUNCS to Prob.ODE
% 050413  joho Added input parameter nEq
% 050418  joho Removed TAB
% 050418  hkh  Removed nEq, is length(Y0), cleaned up  
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050502  hkh  Revised, use input Prob, if given, define full Prob.ODE
% 050705  med  Help updated