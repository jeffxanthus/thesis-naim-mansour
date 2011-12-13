% fminunc is the TOMLAB equivalent to fminunc in Optimization TB
% The unconstrained solver actually used is selectable, see Prob.Solver.Tomlab.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% fminunc solves unconstrained problems of the OPTTB form:
%
%	min                f(x)
%        x
%
% TOMLAB fminunc solves unconstrained problems of the form:
%
%	min                f(x)
%        x
%                 x_L <=   x    <= x_U
%
% OPTTB fminunc does not handle the simple bounds x_L <= x <= x_U
%
% If using a Tomlab solver, the simple bounds could be supplied
% by using the input structure Prob, in fields Prob.x_L and Prob.x_U
%
% function [x, f_k, ExitFlag, Output, g, H, Result] = fminunc(...
%     Func, x_0, options, Prob, varargin)
%
% INPUT: ( 2 arguments always needed )
%
% Func     Function that computes the function f(x) (maybe g(x),H(x))
% x_0      Starting value for the design variables
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%          GradObj, Hessian, HessPattern, HessUpdate, LineSearchType,
%          MaxFunEvals, MaxIter, DiffMinChange and DiffMaxChange,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%
%          Setting HessUpdate to 'steepdesc' gives a CG-algorithm,
%          the Polak-Ribiere algorithm, if running TOMLAB ucSolve
%          Using Prob.Solver.Alg all algorithms are available
%
% Prob     The TOMLAB problem input structure, or the first extra user
%          input argument
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
% Prob.x_L Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% Prob.x_U Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
%
% Additional input as fields in Prob:
%
% Prob.Solver.Tomlab  Name of the Tomlab solver to use instead of the
%                     Optimization Toolbox solver FMINUNC
%
% Prob.Solver.Tomlab = 'snopt';   selects the Tomlab /SOL SNOPT solver
% Prob.Solver.Tomlab = 'ucSolve'; selects the Tomlab ucSolve solver
%
% Default is to use the solver returned by GetSolver('uc', length(x) > 200,0)
%
% See the online help for ucSolve (help ucSolve) and snopt (help snopt and
% help snoptTL. TL is added to the interface file for all MEX solvers).
%
% Prob.Solver.Alg selects different methods in ucSolve (or sTrustr,conSolve).
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.Solver.Tomlab = 'snopt';
%
% OUTPUT:
%
% x        Optimal design parameters
% f_k      Final function value
% ExitFlag exit condition of fminunc.
%      > 0 fminunc converged to a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 Errors, ExitFlag=-Inform, see the Inform parameter description
%
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.funcCount     Number of function evaluations
%   Output.algorithm     Solver name: Type of algorithm used
%   Output.steplength    Length of last step
%   Output.cgiterations  Number of CG iterations (if used)
%   Output.firstorderopt The first-order optimality (if used)
%
% g        The final gradient
% H        The final Hessian
% Result   The Tomlab result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (Tomlab format)
% Structure Result. Fields used (Also see help for conSolve):
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. lower + upper bounds
%   f_k      Function value at x_k
%   g_k      Gradient at x_k
%   H_k      Hessian matrix at x_k
%   Solver   TOMLAB solver, e.g. snopt, knitro or ucSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 29, 1999.  Last modified Jun 6, 2008.

function [x, f_k, ExitFlag, Output, g, H, Result] = fminunc(...
    Func, x_0, options, Prob, varargin)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

global otxProb alphaV

if nargin < 4, Prob = [];
   if nargin < 3, options = [];
      if nargin < 2 
         error('fminunc requires two input arguments Func and x_0');
end, end, end

% ------------Initialization----------------

if ~isempty(Prob) & ~isstruct(Prob)
   Psave=Prob;
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
   Prob.varargin = [{Psave},varargin];
elseif isempty(Prob)
   if ~isempty(otxProb)
      Prob = otxProb;
   else
      Prob = ProbDef;
   end
   if nargin > 11
      Prob.varargin = [{[]},varargin];
   else
      Prob.varargin = varargin;
   end
else
   if isfield(Prob,'TOMLAB')
      Prob.varargin = varargin;
   else
      % Assume structure is part of extra input
      Psave=Prob;
      if ~isempty(otxProb)
         Prob = otxProb;
      else
         Prob = ProbDef;
      end
      Prob.varargin = [{Psave},varargin];
   end
end

z = version;
M7 = z(1) == '7';

solvType=checkType('uc');
Prob.probType = DefPar(Prob,'probType',solvType);

Prob.N   = DefPar(Prob,'N',length(x_0));
n        = Prob.N;
Prob.mLin    = 0;
Prob.mNonLin = 0;
if n == 0
   error('x_0 is empty')
end

Prob.x_L = DefPar(Prob,'x_L',-inf*ones(n,1));
Prob.x_U = DefPar(Prob,'x_U', inf*ones(n,1));

Prob=ProbCheck(Prob,'fminunc',solvType);

if isempty(x_0), x_0 = zeros(n,1); end

Prob.x_0     = x_0(:);
Prob.OPTTB.x = x_0;
Prob.OPTTB.M7= 0;

Prob=ProbCheck(Prob,'fminunc',solvType);

% Set default options

% TolX
Prob.optParam.eps_x=1E-6;

% TolFun
Prob.optParam.eps_f=1E-6;

% DerivativeCheck - off

% Diagnostics - off
Diagnostic=0;

% Hessian/gradient numerical approximation on
Prob.NumDiff=1; 
% No Automatic Differentiation
Prob.ADObj=0;

% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

if n < 200
   % LargeScale off
   Prob.LargeScale=0;
else
   % LargeScale on
   Prob.LargeScale=1;
end

% MaxFunEvals 100*numberOfVariables
% MaxIter 400
Prob.optParam.MaxIter=max(400,100*n);

% DiffMaxChange 1e-1 DiffMinChange 1e-8
Prob.optParam.DiffInt=1E-6;
%Prob.optParam.DiffGradMaxChange=1E-1;

% 'LineSearchType','quadcubic' default, but in Tomlab we set Cubic as default
Prob.LineParam.LineAlg=1;

% PrecondBandWidth 0    TypicalX ones(numberOfVariables,1) 
% MaxPCGIter       max(1,floor(numberOfVariables/2))        TolPCG 0.1

if ~isempty(options)
   Prob.optParam.MaxIter = DefPar(options,'MaxIter',Prob.optParam.MaxIter);
   Prob.optParam.eps_f = DefPar(options,'TolFun',Prob.optParam.eps_f);
   Prob.optParam.eps_x = DefPar(options,'TolX',Prob.optParam.eps_x);
   Prob.optParam.DiffInt = DefPar(options,'DiffGradMinChange',...
        Prob.optParam.DiffInt);

   if isfield(options,'GradObj')
      if strcmpi('on',options.GradObj)
         Prob.NumDiff=0;
      elseif strcmpi('off',options.GradObj)
         Prob.NumDiff=1;
      end
   end
   if isfield(options,'Hessian')
      if strcmpi('on',options.Hessian)
         Prob.NumDiff=0;
      elseif strcmpi('off',options.Hessian)
         if Prob.NumDiff == 0
            Prob.NumDiff=-1;
         end
      end
   end

   Prob.HessPattern = DefPar(options,'HessPattern',Prob.HessPattern);
   Prob.HessMult = DefPar(options,'HessMult',[]);
   if ~isempty(Prob.HessMult)
      fprintf('\nWARNING! Inefficient to use option HessMult with Tomlab\n\n')
      fprintf('\nUse explicit sparse Hessian instead\n\n')
   end

   if isfield(options,'LargeScale')
      if strcmpi('on',options.LargeScale)
         Prob.LargeScale=1;
      elseif strcmpi('off',options.LargeScale)
         Prob.LargeScale=0;
      end
   end
   if isfield(options,'LineSearchType')
      if strcmpi('cubicpoly',options.LineSearchType)
         Prob.LineParam.LineAlg=1;
      elseif strcmpi('quadcubic',options.LineSearchType)
         Prob.LineParam.LineAlg=0;
      end
   end
   if isfield(options,'Display')
      if strcmpi('iter',options.Display)
         Prob.optParam.IterPrint=1;
         PriLev=1;
      elseif strcmpi('final',options.Display)
         Prob.optParam.IterPrint=0;
         PriLev=1;
      elseif strcmpi('off',options.Display)
         Prob.optParam.IterPrint=0;
         PriLev=0;
      end
   end
   if isfield(options,'Diagnostics')
      if strcmpi('on',options.Diagnostics)
         Diagnostic=1;
      elseif strcmpi('off',options.Diagnostics)
         Diagnostic=0;
      end
   end
else
   Prob.HessMult = [];
end

Prob = mkbound(Prob);

if isa(Func,'cell')
   funArgIn  = 1;
   if length(Func) == 1
      Prob.OPTTB.f=Func{1};
      funArgOut = -1;
   elseif length(Func) == 2
      Prob.OPTTB.f=Func{1};
      Prob.OPTTB.g=Func{2};
      if isempty(Prob.OPTTB.g)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      else
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      end
   else
      Prob.OPTTB.f=Func{1};
      Prob.OPTTB.g=Func{2};
      Prob.OPTTB.H=Func{3};
      if isempty(Prob.OPTTB.g)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      elseif isempty(Prob.OPTTB.H)
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      else
         funArgOut = -3;
         Prob.NumDiff=0;
      end
   end
elseif isa(Func,'function_handle')
   Prob.OPTTB.f=Func;
   if M7
      Prob.OPTTB.M7 = 1;
      funArgIn  = abs(nargin(Func));
      funArgOut = abs(nargout(Func));
      if funArgOut == 1
         Prob.NumDiff = max(1,Prob.NumDiff);
      elseif funArgOut == 2
         Prob.NumDiff = min(-1,Prob.NumDiff);
      else
         Prob.NumDiff=0;
      end
   else
      SS = functions(Func);
      funArgIn  = abs(nargin(SS.function));
      funArgOut = abs(nargout(SS.function));
   end
else
   Prob.OPTTB.f=Func;
   if strcmpi(Func,'nlp2_fgH')
      % Now Tomlab calls Tomlab through intermediate fmincon interface - weird
      % Must set proper NumDiff, in case derivatives not defined
      if isempty(Prob.FUNCSX.g)
         funArgOut = 1;
      elseif isempty(Prob.FUNCSX.H)
         funArgOut = 2;
      else
         funArgOut = 3;
      end
      funArgIn  = 3;
   else
      funArgOut = abs(nargout(Func));
      funArgIn  = abs(nargin(Func));
   end
   if funArgOut == 1
      Prob.NumDiff = max(1,Prob.NumDiff);
   elseif funArgOut == 2
      Prob.NumDiff = min(-1,Prob.NumDiff);
   else
      Prob.NumDiff=0;
   end
end

Prob.OPTTB.funArgOut=funArgOut;
Prob.OPTTB.funArgIn =funArgIn;

% Must save r and J if the underlying problem is a least squares problem
r = Prob.FUNCS.r;
J = Prob.FUNCS.J;

Prob = tomFiles(Prob,'optim_fgH','optim_g','optim_H');

Prob.FUNCS.r=r;
Prob.FUNCS.J=J;

Prob.Solver.Alg  = DefPar(Prob.Solver,'Alg',0);
Solver           = DefPar(Prob.Solver,'Tomlab',[]);

if isempty(Solver)
   Solver = GetSolver('uc',n > 200 | Prob.LargeScale,0);
   Prob.Solver.Tomlab = Solver;
else
   fprintf('TOMLAB Solver: %s\n',Solver);
end

if Diagnostic
   Alg=Prob.Solver.Alg;
   if isempty(Alg)
      Alg=1;
   end
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf('   Diagnostic Information');
   fprintf('\n');
   if ~isempty(Prob.Name)
      fprintf('Problem %s',Prob.Name);
   end
   fprintf('\n');
   fprintf('Number of variables: %d\n',n);
   fprintf('Functions\n');
   fprintf(' Objective: TOMLAB optim_fgH is calling  %s',Func);
   fprintf('\n');
   fprintf(' Gradient:  TOMLAB optim_g is');
   if Prob.NumDiff==0
      fprintf(' using gradient in global variables NLP_g\n');
   else
      fprintf(' using finite-differencing method # %d\n',Prob.NumDiff);
   end
   fprintf('\n');
   if Prob.NumDiff==0
      fprintf(' Hessian:   TOMLAB optim_H is');
      fprintf(' using Hessian in global variables NLP_H\n');
   elseif Prob.NumDiff>0 & strcmpi(Solver,'ucSolve') & Alg==0
      fprintf(' Hessian:   TOMLAB optim_H is');
      fprintf(' using finite-differencing method # %d\n',Prob.NumDiff);
   elseif strcmpi(Solver,'ucSolve') & Alg > 0
      fprintf(' Hessian:   TOMLAB is');
      fprintf(' using quasi-Newton update\n');
   elseif strcmpi(Solver,'sTrustr')
      fprintf(' Hessian:   TOMLAB is');
      fprintf(' using quasi-Newton update\n');
   end
   fprintf('\n');
   fprintf(' Number of lower bound constraints:              %d\n',...
           sum(~isinf(Prob.x_L)));
   fprintf(' Number of upper bound constraints:              %d\n',...
           sum(~isinf(Prob.x_U)));
   fprintf('\n');
   fprintf('Solver: %s',Solver);
   if isempty(Prob.Solver.Alg)
      fprintf('. Default algorithm.\n');
   else
      fprintf('. Alg # %d.\n',Prob.Solver.Alg);
   end
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf(' End diagnostic information\n');
   fprintf('\n');
end

if strcmpi('sTrustr',Solver)
   Result   = tomRun('sTrustr',Prob);
elseif strcmpi('ucSolve',Solver)
   if isfield(options,'HessUpdate')
      if strcmpi('bfgs',options.HessUpdate)
         Prob.Solver.Alg=1;
      elseif strcmpi('dfp',options.HessUpdate)
         Prob.Solver.Alg=4;
      elseif strcmpi('steepdesc',options.HessUpdate)
         Prob.Solver.Alg=6;
      end
   elseif isempty(Prob.Solver.Alg)
      Prob.Solver.Alg=1; % Default BFGS
   end
   Result   = tomRun('ucSolve',Prob);
elseif strcmpi('sTrustr',Solver)
   Result   = tomRun('conSolve',Prob);
else
   Result   = tomRun(Solver,Prob);
end

PrintResult(Result,Prob.PriLev);

x        = x_0;          % Initialize x to x_0 to get correct size
x(:)     = Result.x_k;   % Column vector Result.x_k gets into x with size(x_0)
f_k      = Result.f_k;
g        = Result.g_k;
H        = Result.H_k;
ExitFlag = Result.ExitFlag;
Inform   = Result.Inform;

% Convert ExitFlag to OPTIM TB
if ExitFlag==0
   ExitFlag=1;
   if PriLev > 0
      fprintf('fminunc (%s',Solver);
      fprintf('): Convergence\n');
   end
else
   switch Inform
      case 101 
        ExitFlag=0;
        if PriLev > 0
           fprintf('fminunc (%s',Solver);
           fprintf('): Too many iterations\n');
        end
      otherwise
        ExitFlag=-Inform;
        if PriLev > 0
           fprintf('fminunc (%s',Solver);
           fprintf('): Possible Error.');
           ExitText = Result.ExitText;
           if ~isempty(ExitText)
              fprintf(' Solver information for Inform = %d:\n',Inform);
              for i = 1:size(ExitText,1)
                  fprintf('%s\n',ExitText(i,:));
              end
           else
              fprintf(' See solver help: Inform = %d\n',Inform);
           end
           fprintf('\n');
        end
   end
end

if nargout > 3
   Output.iterations  = Result.Iter;
   Output.funcCount   = Result.FuncEv;
   Output.algorithm   = [ Solver ': ' Result.SolverAlgorithm];
   if isempty(alphaV)
      Output.stepsize   = 1;
   else
      Output.stepsize   = alphaV(length(alphaV));
   end
   Output.firstorderopt= [];
end

% MODIFICATION LOG
%
% 011203 hkh General solver handling with Prob.Solver.TOM
% 011204 hkh Test is nonempty before setting value from options into Prob
% 011204 hkh Save r and J for least squares problems
% 030126 hkh Improve handling of lower and upper bounds
% 030127 hkh Handle case when x matrix
% 030128 hkh Use DefPar, improve comments, use Prob.N to mark a structure
% 031201 hkh Change AutoDiff to field ADObj
% 040203 hkh Use field FuncEv instead of global n_f
% 040203 hkh Print ExitText if possible error (and PriLev > 0)
% 050422 hkh Handle new type of function handle in Matlab 7.x
% 050422 hkh Test and set NumDiff if calling Tomlab solver
% 060217 hkh Set LargeScale=1 as default only if Prob.N >= 200
% 060217 hkh Lower Prob.DiffInt from 1E-8 to 1E-6
% 060219 hkh Print TOMLAB Solver if the solver is changed from default
% 060814 med FUNCS used for callbacks instead
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added