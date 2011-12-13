% plotInit.m:
%
% Initialization routine that sets up the pParam structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Oct 11, 1998.    Last modified Nov 6, 2000.

function pParam = plotInit(Prob, pParam)

if nargin < 2
   pParam = [];
end

pParam.Prob=Prob;

pParam.ask = -1;
pParam.PriLev = Prob.PriLevOpt; % Print Level in driver routine
pParam.File = 'test';      % generate file name as a string (?.m and ?.mat): 

pParam.plotType = 1;       % Type of plot
pParam.viewPnt= [-37.5 30];% View point for mesh plot
pParam.ixAxis = [1 2];     % Subspace to plot if dimension > 2 
pParam.PlotPnts=[60 60];   % Number of points used in twodimensional grid

pParam.HoldRun = 0;        % =1 if "hold previous run"

n = Prob.N;

if isempty(n) | n == 0
   n = max([length(Prob.x_0), length(Prob.x_L), length(Prob.x_U)]);
   Prob.N = n;
end
% Axis for contour and mesh plot. 

if isfield(Prob,'x_max')
   x_max = Prob.x_max;
else
   x_max = [];
end

if isempty(x_max) 
   if isfield(Prob,'x_U')
      x_U = Prob.x_U;
      if isempty(x_U), x_U = Inf*ones(n,1); end
      if isinf(x_U(1)) 
         if isempty(Prob.x_0)
            x_max(1) = 10;
         else
            x_max(1) = max(5,Prob.x_0(1)+1);
         end
      else
            x_max(1) = x_U(1);
      end
      if n > 1
         if isinf(x_U(2)) 
            if isempty(Prob.x_0)
               x_max(2) = 10;
               else
               x_max(2) = max(5,Prob.x_0(2)+1);
            end
         else
               x_max(2) = x_U(2);
         end
      end
   else
      x_max = [10;10];
   end
end

if isfield(Prob,'x_min')
   x_min = Prob.x_min;
else
   x_min = [];
end

if isempty(x_min)
   if isfield(Prob,'x_L')
      x_L = Prob.x_L;
      if isempty(x_L), x_L = -Inf*ones(n,1); end
      if isinf(x_L(1)) 
         if isempty(Prob.x_0)
            x_min(1) = -10;
         else
               x_min(1) = min(-5,Prob.x_0(1)-1);
         end
      else
         x_min(1) = x_L(1);
      end
      if n > 1
         if isinf(x_L(2)) 
            if isempty(Prob.x_0)
               x_min(2) = -10;
               else
               x_min(2) = min(-5,Prob.x_0(2)-1);
            end
         else
            x_min(2) = x_L(2);
         end
      end
   else
      x_min = [-10;-10];
   end
end
n1 = length(x_min);
n2 = length(x_max);
x_min = x_min(:);
x_max = x_max(:);
if n1 ~= n2
   if n1 > n2
      x_max = [x_max;x_min(n2+1:n1)];
   else
      x_min = [x_min;x_max(n1+1:n2)];
   end
end

ix = x_min == x_max;
x_min(ix)=x_min(ix)-1;

pParam.xMin1 = x_min(1);     
pParam.xMax1 = x_max(1);     
if n > 1
   pParam.xMin2 = x_min(2);    
   pParam.xMax2 = x_max(2);   
else
   pParam.xMin2 = [];
   pParam.xMax2 = [];
end
if nargin < 2
   pParam.beta1   = [];
   pParam.beta2   = [];
   pParam.p_dx    = [];
   pParam.alphaV  = [];
   pParam.F_X     = [];
end

pParam.Solver     = Prob.Solver.Name;
pParam.Alg        = Prob.Solver.Alg;
pParam.Method     = Prob.Solver.Method;
pParam.meshXcenter= Prob.x_0;
pParam.N          = length(Prob.x_0);
pParam.X_min      = Prob.x_min;
pParam.X_max      = Prob.x_max;
%pParam.X_min     = x_min;
%pParam.X_max     = x_max;
pParam.figs       = 2;

% MODIFICATION LOG:
%
% 981011 hkh  Defined routine from code in xxxOpt3.m. Changed fields,
%             including Prob in pParam. Also setting all other fields
%             to default values
% 981022 hkh  Add definition of pParam.N
% 981027 hkh  Add pParam.Solver and get alg from Prob.Solver.Alg
% 981110 hkh  Eliminate MAX_x,MAX_c,MAX_r from pParam
% 000927 hkh  Major revison, bugs in handling of upper and lower