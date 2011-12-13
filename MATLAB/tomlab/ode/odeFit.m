% odeFit - General driver routine for ode fitting in TOMLAB
%
% If using the TOMLAB format (TQ), call with:
%
%   function [Result] = odeFit(Solver, odeSolver, Prob, PriLev, ask);
%
% The following call will also work (similar 7-input format as below)
%
%   function [Result] = odeFit(Solver, odeSolver, [], [], Prob, PriLev, ask);
%
% If using the TOMLAB Init File format, call with:
%
%   function [Result] = odeFit(Solver, odeSolver, probFile, probNumber, ...
%                              Prob, PriLev, ask);
%
% A third alternative is the call
%
%   function [Result] = odeFit(Solver, odeSolver, probNumber, Prob, PriLev,
%   ask);
%
% Then the default file for ode fitting problems is used.
%
% If calling with odeFit; (no arguments), a list of available solvers is given
%
% If calling with odeFit(Solver);
%    a list of available solvers for solving odes is given
%
% tomRun checks if the third argument is a string, a structure, a number,
% or is empty, to determine which input format is used.
%
% if isempty(Solver),    odeFit is using the TOMLAB default solver for
%                        ode fitting as default
% if isempty(odeSolver), odeFit is using the TOMLAB default odesolver for
%                        ode fitting as default
% if isempty(probFile),  no default probFile availible yet
%
%
% A problem available in the TOMLAB Init File format is defined using a call
%           Prob=probInit(probFile, probNumber, ask, Prob)
%
% INPUT: (if [] is given or less parameters are given, default values are used)
%
% Solver     The name of the solver that should be used to optimize the problem.
%            If the Solver may run several different optimization algorithms,
%            then the values of Prob.Solver.Alg and Prob.Solver.SubAlg
%            determines which algorithm.
%            The Solver name is put in Prob.Solver.Name
%
% odeSolver The name of the ordinary differential equations solver used to
%            optimize the problem. If the Solver may run several different
%            optimization algorithms, then the values of Prob.?ODE.Alg
%            and Prob.?ODE.SubAlg determines which algorithm.
%
% probFile   User problem initialization file.
%
% probNumber Problem number in probFile.
%            If empty of left out, either probNumber=Prob.P (if set) or
%            otherwise probNumber=1.
%            When calling the probFile with probNumber=0, probFile must
%            return a string matrix with the names of the problems defined.
%
% Prob       Problem structure. Either define the structure with the
%            call:  Prob=probInit(probFile, probNumber, ask);
%            or set in Prob the parameters with special values.
%            See the manual for a description of the Prob structure
%
%            P        probNumber (YOU MUST SET THIS, IF SETTING PROB AS INPUT!)
%
%            Examples of other fields to set:
%            probFile probFile
%            uP       User problem parameters, which you can use when computing
%                     the functions. See the variable ask.
%            optParam A substructure with optimization parameters. Default value
%                     from optParamSet(Solver,probType)
%            x_0      Starting point for the problem.
%
%            ODE      Struct with the ordinary differential equations
%                     subproblems.
%
%                     E         Matrix with the measured data.
%                     t         Vector with the times when the
%                               data was measured.
%                     resWeight Vector with weights on the
%                               residuals.
%                     weight    the weight of each data series
%
% PriLev     Print level when displaying the result of the optimization in
%            routine PrintResult.
%            =0 No output, =1 Final result, shorter version,
%            =2 Final result, longer version, =3 All output
%            If isempty(PriLev), Prob.PriLev is used, if nonempty
%
%            The printing level in the optimization solver is controlled
%            by setting the parameter Prob.PriLevOpt
%
% ask        If ask>=1: ask questions in probFile. If ask = 0: use defaults
%            If ask < 0: use values in user params uP if defined or defaults.
%            If isempty(ask): If length(uP) > 0, ask=-1, else ask=0
%
% OUTPUT:
% Result Structure with optimization results. See manual for description
%        Result.Prob holds the input Prob structure

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Apr 21, 2005.  Last modified Oct 30, 2008.

function [Result] = odeFit(Solver,odeSolver, Prob, P4, P5, P6, P7)

if nargin < 3    
   if nargin < 1
      for i=1:10
          fprintf('\nSolvers for problem type %d\n',i)
          z1=SolverList(i);
          disp(z1);
      end
      return
   elseif nargin == 1
      if ~ischar(Solver)
         fprintf('\nSolvers for problem type %d\n',Solver)
         z1=SolverList(max(1,min(10,round(Solver))));
         disp(z1);
      else
         for i=1:10
             fprintf('\nSolvers for problem type %d\n',i)
             z1=SolverList(i);
             disp(z1);
         end         
      end
      for i=1:1
         fprintf('\nSolvers for ode-problems:\n')
         y=odeSolverList(i);
         disp(y);
      end
      return
   else
      fprintf('\n');
      error('odeFit needs at least three parameters');
   end    
end

if nargin < 7
   P7=[];
   if nargin < 6
      P6=[];
      if nargin < 5
         P5=[];
         if nargin < 4
            P4=[];
end, end, end, end

if isempty(Prob)
   error('Prob must not be empty');
elseif ~isstruct(Prob)
   error('Prob must be of type struct'),
end

% Set default odeSolver
if isempty(odeSolver)
   odeSolver = 'lsode';
else
   odeSolver = deblank(odeSolver);
end

if isempty(Solver)
   if strcmpi(odeSolver,'modfit')
       Solver = 'nlpql';
   else
       Solver = 'SNOPT';
   end
else
   Solver=deblank(Solver);
end

% Check modfit solver combinations.
switch lower(Solver)
  case {'dfnlp','dn2gb','dslmdf','dfneld'}
    modSolver = 1;
  otherwise
    modSolver = 0;
end
switch lower(odeSolver)
  case {'dopri5','radau5','ind-dir'}
     % OK!
  otherwise
    if modSolver == 1
       fprintf('Odesolver %s NOT allowed.',odeSolver);
       fprintf(' For a modfit solver, odesolver must be any of:\n')
       fprintf(' dopri5, radau5 or ind-dir\n\n')
       error('Illegal solver choice!')
    end
end

switch lower(odeSolver)
    case {'lsode','rksuite',...
         'ode23','ode113','ode15s','ode23s','ode23t','ode23tb','ode45',...
         'dopri5','radau5','ind-dir'}
         Prob.SolverODE = odeSolver;
    otherwise
         fprintf('odeSolver %s',odeSolver);
         fprintf(' NOT found\n')
         error('Illegal odeSolver choice!')
end

% Check if solver is OK

switch lower(Solver)
    case {'snopt','clssolve','npsol','nlssol','minos','nlpsolve',...
         'filtersqp','fmincon','consolve','snopt7','snopt6','nlpql',...
         'dfnlp','dn2gb','dslmdf','dfneld',...
         'oqnlp','knitro','lsgrg2','conopt','lgo'}
    otherwise
         fprintf('Solver %s',Solver);
         fprintf(' is not possible to use for ODE parameter estimation\n')
         error('Illegal ODE Optimization Solver!')
end

Prob.probType = checkType('ode');

if modSolver 
   % Run modfit 
   % Check ODE solution specific fields in Prob.ODE
   Prob             = odeProbCheck(Prob);
   Prob.optParam    = optParamDef(Solver,checkType('cls'),length(Prob.x_0),0,0);
   Prob.Solver.Name = Solver;
   Result = modfitTL(Prob);
   if P4 > 0
      PrintResult(Result,P4);
   end
else
   Result = tomRun(Solver, Prob, P4, P5, P6, P7);
end

% MODIFICATION LOG:
%
% 050412  bkh  Created from copy of tomRun.m
% 050413  bkh  First untested version ready
% 050417  hkh  Revised, pure interface calling tomRun, check Solver is OK
% 050418  hkh  Add call to modfitTL
% 050419  bkh  Add call to all ML odesolvers except ode15i
% 050429  bkh  Split modfit into separate solvers and odesolvers
% 050429  bkh  Enabled use of modfit odesolvers with any solver
% 050503  hkh  Call odeProbCheck before modfitTL
% 050602  hkh  Use checkType to get probType for ode
% 050705  med  Help updated
% 050801  med  isstr replaced by ischar
% 081031  ango snopt7 replaces snopt6