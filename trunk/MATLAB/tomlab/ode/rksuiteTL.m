% TOMLAB RKSUITE ODE Solver
%
% RKSUITE is a suite of codes based on Runge-Kutta formulas that solves the
% initial value problem for a first order system of ordinary differential
% equations using Runge-Kutta methods.
%
% Some about catastrofic failures are redirected from stdout to the
% file out.txt. 
%
%  Fields used in struct Prob
%  
%   ODE       Structure containing data specific to ODE:s   
%   
%   PriLev    For security reasons, a file out.txt is always opened so that
%             rksuite has somewhere to place outputs. 
% 
%             PriLev <= 0 Only print failures that causes RKSUITE to
%             stop.          
%             PriLev  = 1 Also print warnings and a short introduction
%             PriLev >= 2 Print all messages
%             Default value is PriLev = 0
%
%   FUNCS.f   Name of m-file that contains the system y' = f(t,y) with
%             t, y, Prob as order of inputs.
%		
%
%  Fields used in struct Prob.ODE  
%
%   Y0        Vector with the initial values of each ODE variable y.
%             nEq = length(Y0) is the number of equations in the system of ODEs.
% 
%   tInit     The initial value of the independent variable.
%  
%   tStop     The integration proceeds from tInit in the direction
%             of tStop. You can, and often will, terminate the integration
%             before reaching tStop, but you cannot integrate past tStop.
%           
%             Constraint: tStop must be clearly distinguishable from tInit
%             in the precision available.
%
%   tWant     Array of points where rksuite will compute solutions. 
%
%   InitStep  InitStep = 0    The code will select automatically the first
%                         step size (recommended choice).
%             InitStep != 0  non-zero  - The code will try ABS(INITSTEP) for
%                         the first step size of the integration.
%       
%             Default value is InitStep = 0
%
%   absTol    Scalar or array of of the same length
%             as the number of equations in the system. 
%             If absTol is scalar, the value of absTol is taken for
%             all absTol(L).     
%             Choose absTol so that the value of a solution
%             component Y(L) is not important when Y(L) is smaller in
%             magnitude than absTol/relTol. 
%             absTol(L) must be greater or equal than relTol times
%             the squareroot of the smallest possible number on the
%             machine being used.  
%	   
%             Default value is absTol(L) = 1e-14
%
%
%
%   relTol    Defines the relative error tolerance.
%             Constraint: 0.01 >= TOL >= 10 * (machine epsilon) 
%             Default relTol value is 1e-7
%
%
%   RKSUITE   Structure containing rksuite specific settings   
%
%
%  Fields used in input structure Prob.ODE.RKSUITE
%
%   
%   METHOD    Specifies which Runge-Kutta formula pair is to be used for
%	          integration. 
%             METHOD = 1 - use the (2,3) pair
%             METHOD = 2 - use the (4,5) pair
%             METHOD = 3 - use the (7,8) pair
%             Constraint: 1 <= METHOD <= 3
%	          Default method is 1
%
%   
%   TASK      Defines which mode rksuite will use.
%             TASK = 0, use UT (Usual Task)
%
%             TASK != 0, use CT (Complicated task)          
%             Step from tInit towards tStop, accepting answers at the
%             points chosen by rksuite. This is often the best way to
%             proceed when you want to see how the solution behaves
%             throughout the interval of integration because the code
%             will tend to produce answers more frequently where the
%             solution  changes more rapidly (the step sizes tend to be
%             smaller there). 
%             TWANT is not concidered in CT, and may thus be left empty.  
%                 
%             The default value is TASK = 0, use UT.
%
%              
%
%   ERRASS    =  0 - do not attempt to assess the true error.
%             =! 0 - assess the true error, the difference between
%                  the numerical solution and the true solution.
%                  (The cost of this is roughly twice the cost 
%                   of the integration itself with METHODs 2 and 3,
%                   and three times with METHOD = 1.)
%             Default value is 0.
%
%
%  Fields used in output structure Result
%
%   ODE       Structure containing ODE specific result.
%
%   Inform    Some kind of Global state of RKSUITE
%
%   ExitFlag  Tomlab defined exit flag
%
%   ExitText  Status text from RKSUITE
%
%   
%
%  Fields defined in structure Result.ODE
%
%  
%   y         Matrix where each column is the solution
%             to the t in corresponding column of T.
%
%   t         Points where rksuite actually computed a solution to the
%             system. 
% 
%   yp        Matrix where each column is the first
%             derivative of the corresponding solution in y.
%   
%  Fields in Result.ODE.RKSUITE. containing RKSUITE specific results.
%
%   yMax      Array with the same number of elements as the number of
%             equations in the system.
%             yMax(L) is the largest value of abs(y(L)) computed at 
%             any step in the integration so far. (With METHODs 1 and 2,
%             y(L) is computed by interpolation, so yMax(L) might be 
%             a little larger in magnitude than any value abs(y(L))
%             reported so far.) 
%
%   Inform    For each point in Result.ODE.t, this is the status given
%             by RKSUITE.  

% These are not returned in RKSUITE NOW ????
%   FuncEv    Number of function evaluations. 
%
%   JacEv     Number of Jacobian evaluations.
%

% Johan Holmgren, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Apr 7, 2005. Last modified Aug 14, 2006.

function Result = rksuiteTL(Prob)

if nargin < 1 
   error('rksuiteTL needs the Prob structure as input');
end

RKSUITE  = DefPar(Prob.ODE, 'RKSUITE', []);
tStop    = DefPar(Prob.ODE, 'tStop', inf);
InitStep = DefPar(Prob.ODE, 'InitStep', 0);

nEq            = length(Prob.ODE.Y0);
RKSUITE.relTol = DefPar(Prob.ODE, 'relTol', 1e-4);
absTol         = DefPar(Prob.ODE, 'absTol', 1e-8);
RKSUITE.absTol = absTol*ones(1,nEq);
Result.Solver  = 'RKSUITE';

% ERROR IN MEX; must put in this for now
Prob.ODE.nEq   = nEq;

[Result.ODE.y, Inform, Result.ODE.t, Result.ODE.yp, Result.ODE.RKSUITE.yMax,...
        Result.ODE.RKSUITE.Inform] = rksuite(Prob.ODE.Y0, ...      
        Prob.ODE.tInit, tStop, Prob.ODE.tWant, InitStep, ...
        Prob.ODE.f, Prob.PriLevOpt, RKSUITE, Prob);

Result.Inform = Inform;

switch Inform
  case 1
   Result.ExitFlag = 0;   
   Result.ExitText = 'Complete success';
  case 2
   Result.ExitFlag = 0;
   Result.ExitText = 'Success';
  case 3
   Result.ExitFlag = 0;
   Result.ExitText = 'Success';
  case 4
   Result.ExitFlag = 0;
   Result.ExitText = 'Problem is stiff';
  case 5
   Result.ExitFlag = 1;
   Result.ExitText = 'Cannot reach the requested accuracy';
  case 6
   Result.ExitFlag = 1;
   Result.ExitText = 'global error not reliable';
  otherwise
   Result.ExitFlag = 1;
   Result.ExitText = 'Catastrofic failure';
end
  
if Result.Inform < 5
   Result.ExitFlag = 0;
else
   Result.ExitFlag = 1;
end

% MODIFICATION LOG:
%
% 050407  joho Written
% 050412  joho Added documentation
% 050418  joho Removed TABS
% 050418  hkh  Changed fields, remove Result.Prob setting, cleaned up
% 050419  joho Changed call from rksuiteMex to rksuite
% 050419  joho Minor changes of documentation
% 050420  joho Moved output Result.Inform 
% 050422  bkh  Changed odeH_s/hStart, odeT_s/tStart and odeT_e/tEnd to 
%              InitStep, tInit and tStop respectively
% 050502  hkh  Clean up, use RKSUITE, not Prob.ODE.RKSUITE
% 060814  med  FUNCS used for callbacks instead