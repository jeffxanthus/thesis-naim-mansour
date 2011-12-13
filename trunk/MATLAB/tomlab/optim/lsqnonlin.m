% lsqnonlin is the TOMLAB equivalent to LSQNONLIN in Optimization TB 
% The NLLS/NLP solver actually used is selectable, see Prob.SolverQP.
% If no active choice of solver is made, GetSolver will be called to select 
% the best solver depending on the license used
%
% Opt tlbx lsqnonlin solves nonlinear least squares problems on the form:
%
%	min       0.5 * sum {r(x).^2} <==>      min       0.5 *  r(x)^T * r(x)
%    x                                       x
%             x_L <=   x  <= x_U                      x_L <=   x  <= x_U
%
% TOMLAB lsqnonlin solves nonlinear least squares problems on the form:
%
%	min       0.5 * sum {r(x).^2} <==>      min       0.5 *  r(x)^T * r(x)
%    x                                       x
%             x_L <=   x  <= x_U                      x_L <=   x  <= x_U
%             b_L <=  Ax  <= b_U                      b_L <=  Ax  <= b_U
%
% where r is the m-dimensional residual vector and
% x is a n-dimensional unknown parameter vector, with bounds x_L and x_U.
%
% The linear constraints are treated by the Tomlab solvers as standard
% but they are not treated in Optimization TB. 
% They could be added in the Tomlab Prob structure, input as extra
% argument or as global variable otxProb, see below
%
% function [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] = lsqnonlin(...
%           rFunc, x_0, x_L, x_U,  options, Prob, varargin)
%
% INPUT: ( 2 arguments always needed )
%
% rFunc    Function that computes the residual r(x) (and maybe J(x))
% x_0      Starting value for the design variables
% x_L      Lower bounds on the design values. -Inf == unbounded below.
%          Empty x_L ==> -Inf on all variables
% x_U      Upper bounds on the design values.  Inf == unbounded above.
%          Empty x_U ==>  Inf on all variables
% options  Replaces the default optimization parameters
%          Fields used: Display, TolX, TolFun, DerivativeCheck, Diagnostics,
%          Jacobian, JacobPattern, LineSearchType, LevenbergMarquardt,
%          MaxFunEvals, MaxIter, DiffMinChange and DiffMaxChange,
%          LargeScale, MaxPCGIter, PrecondBandWidth, TolPCG, TypicalX.
%          JacobMult: If nonempty, Tomlab makes i=1:dim(x) calls
%             feval(JacobMult,J,e_i,-1,varargin),   to obtain the J matrix.
%          e_i= (....,1, ....), 1 in the ith position, 0:s otherwise
%          The Tomlab solvers cannot utilize the J*x computation in JacobMult
%
%    NOTE: If Jacobian is on, the user defines the Jacobian matrix
%          and returns it together with the residual vector in rFunc
%          [r_k,J_k]=feval(rFunc,x,Prob)
%          The Prob extra argument is possible in the TOMLAB lsqnonlin,
%          and can be used to send information to rFunc.
%          If r_k is in R^m and x in R^n, then J_k is an m-by-n matrix
%          where J(i,j) is the partial derivative of r_k(i) with respect to x(j)
%
% Prob     The TOMLAB problem input structure, or the first extra user
%          input argument
%          If defining your own limited Tomlab input structure, first do
%             Prob = ProbDef;
%          Then set fields in this structure
%
% Additional fields used in the Prob structure, beside standard fields:
%
% Prob.Solver.Tomlab  Name of the Tomlab solver to use instead of the
%                     Optimization Toolbox solver LSQNONLIN
%
% Prob.Solver.Tomlab ='nlssol';   selects the TOMLAB /SOL nlssol solver
% Prob.Solver.Tomlab ='clsSolve'; selects the TOMLAB clsSolve solver
% Prob.Solver.Tomlab ='snopt';    selects the TOMLAB /SOL snopt§ solver
% Prob.Solver.Alg 1-4 selects the 4 different methods in clsSolve.
%
% Default is to use the solver returned by:
%            GetSolver('cls', length(x) > 200 | Prob.LargeScale,0)
%
% See the online help for cls§olve (help cls§olve) and nlssol (help nlssol and
% help nlssolTL. TL is added to name of interface file for all MEX solvers).
%
% Note!  Another way to input the Tomlab problem structure is to define
% a global structure, called otxProb, with any of the fields in the Tomlab
% Prob input structure format that you want to set. Do
%          global otxProb
%          otxProb = ProbDef;  % Create an empty Tomlab Prob structure
%          "set fields in otxProb", e.g. otxProb.Solver.Tomlab = 'nlssol';
%
% OUTPUT:
%
% x        Optimal design parameters
% f_k      Optimal residual sum of squares, sum {r(x).^2} (Note! no 0.5)
% r_k      The optimal residual vector
% ExitFlag exit condition of lsqnonlin.
%      > 0 lsqnonlin converged to a solution X.
%        0 Reached the maximum number of iterations without convergence
%      < 0 Errors, ExitFlag=-Inform, see the Inform parameter description
%
% Output   Structure. Fields:
%   Output.iterations    Number of iterations
%   Output.funcCount     Number of function evaluations
%   Output.algorithm     Type of algorithm used
%   Output.steplength    Length of last step
%   Output.cgiterations  Number of CG iterations (if used)
%   Output.firstorderopt The first-order optimality (if used)
% Lambda   Structure with Lagrange multipliers at the solution
%    Lambda.lower     Lagrange multipliers for the lower bounds lb
%    Lambda.upper     Lagrange multipliers for the upper bounds ub
% J_k      The optimal Jacobian matrix
% Result   The TOMLAB result output structure
%
% ADDITIONAL OUTPUT PARAMETERS (TOMLAB format)
% Structure Result. Fields used (Also see help for clsSolve):
%   Iter     Number of iterations
%   ExitFlag Exit flag
%            == 0  => OK
%   Inform   If ExitFlag > 0, Inform=ExitFlag.
%   x_k      Solution
%   v_k      Lagrange parameters. Constraints + lower + upper bounds
%   f_k      Function value 0.5*x'*F*x+c'*x
%   r_k      Residual at x_k
%   g_k      Gradient at x_k
%   J_k      Jacobian matrix at x_k
%   Solver   TOMLAB solver, e.g. nlssol, snopt or clsSolve
%   SolverAlgorithm  Description of method used
%   x_0      Starting point x_0
%   xState   State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 2, 1999.   Last modified Jun 6, 2008.

function [x, f_k, r_k, ExitFlag, Output, Lambda, J_k, Result] = lsqnonlin(...
  rFunc, x_0, x_L, x_U,  options, Prob, varargin)

if nargin == 1 & strcmpi(f,'defaults')
    x = struct;
    return
end

global alphaV n_r otxProb

if nargin < 6, Prob = [];
   if nargin < 5, options = [];
      if nargin < 4, x_U = []; 
         if nargin < 3, x_L = []; 
            if nargin < 2 
	       error('lsqnonlin requires two input arguments rFunc and x');
end, end, end, end, end

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

n = max([length(x_0(:)),length(x_L(:)),length(x_U(:))]);

Prob.OPTTB.M7= 0;
Prob.OPTTB.x = x_0;
Prob.x_0     = x_0(:);
Prob.x_L     = x_L(:);
Prob.x_U     = x_U(:);
Prob.N = n;
Prob.mLin = 0;
Prob.mNonLin = 0;

solvType=checkType('cls');
Prob.probType=checkType('ls');

Prob = tomFiles(Prob,'ls_f','ls_g','ls_H',[],[],[],'optim_rJ','optim_J');

Prob=ProbCheck(Prob,'lsqnonlin',solvType);

if ~isempty(Prob.A)
   Prob.probType=solvType;
end

Prob = mkbound(Prob);

% Set default options

% TolX
Prob.optParam.eps_x=1E-6;

% TolFun
Prob.optParam.eps_f=1E-6;

% DerivativeCheck - off

% Diagnostics - off

Diagnostic=0;

% Jacobian off
Prob.NumDiff=1; 
% No Automatic Differentiation
Prob.ADObj=0;

% Display final is default
Prob.optParam.IterPrint=0;
PriLev=1;

% LargeScale on
Prob.LargeScale=1;

% MaxFunEvals 100*numberOfVariables

% MaxIter 400
Prob.optParam.MaxIter=max(400,100*n);

% DiffMaxChange 1e-1 DiffMinChange 1e-8
Prob.optParam.DiffInt=1E-8;
%Prob.optParam.DiffGradMaxChange=1E-1;

% 'LineSearchType','quadcubic' default, but in TOMLAB we set Cubic as default
Prob.LineParam.LineAlg=1;

% LevenbergMarquardt    on

% PrecondBandWidth 0    TypicalX ones(numberOfVariables,1) 
% MaxPCGIter       max(1,floor(numberOfVariables/2))        TolPCG 0.1

if ~isempty(options)
   Prob.optParam.MaxIter = DefPar(options,'MaxIter',Prob.optParam.MaxIter);
   Prob.optParam.eps_f = DefPar(options,'TolFun',Prob.optParam.eps_f);
   Prob.optParam.eps_x = DefPar(options,'TolX',Prob.optParam.eps_x);
   Prob.optParam.DiffInt = DefPar(options,'DiffGradMinChange',...
        Prob.optParam.DiffInt);

   if isfield(options,'Jacobian')
      if strcmpi('on',options.Jacobian)
         Prob.NumDiff=0;
      elseif strcmpi('off',options.Jacobian)
         Prob.NumDiff=1;
      end
   end
   if isfield(options,'JacobPattern')
      if ~isempty(options.JacobPattern)
         Prob.JacPattern=options.JacobPattern;
      end
   end
   Prob.JacobMult = DefPar(options,'JacobMult',[]);
   if ~isempty(Prob.JacobMult)
      fprintf('\nWARNING! Inefficient to use option JacobMult with Tomlab\n\n')
      fprintf('\nUse explicit sparse Jacobian instead\n\n')
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
   Prob.JacobMult = [];
end

Prob = mkbound(Prob);

if isa(rFunc,'cell')
   funArgIn  = 1;
   if length(rFunc) == 1
      Prob.OPTTB.r=rFunc{1};
      funArgOut = -1;
   elseif length(rFunc) == 2
      Prob.OPTTB.r=rFunc{1};
      Prob.OPTTB.J=rFunc{2};
      if isempty(Prob.OPTTB.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      else
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      end
   else
      Prob.OPTTB.r=rFunc{1};
      Prob.OPTTB.J=rFunc{2};
      Prob.OPTTB.d2r=rFunc{3};
      if isempty(Prob.OPTTB.J)
         Prob.NumDiff = max(1,Prob.NumDiff);
         funArgOut = -1;
      elseif isempty(Prob.OPTTB.d2r)
         Prob.NumDiff = min(-1,Prob.NumDiff);
         funArgOut = -2;
      else
         funArgOut = -3;
         Prob.NumDiff=0;
      end
   end
elseif isa(rFunc,'function_handle')
   Prob.OPTTB.r=rFunc;
   if M7
      Prob.OPTTB.M7 = 1;
      funArgIn  = abs(nargin(rFunc));
      funArgOut = abs(nargout(rFunc));
      if funArgOut == 1
         Prob.NumDiff = max(1,Prob.NumDiff);
      elseif funArgOut == 2
         Prob.NumDiff = min(-1,Prob.NumDiff);
      else
         Prob.NumDiff=0;
      end
   else
      SS = functions(rFunc);
      funArgIn  = abs(nargin(SS.function));
      funArgOut = abs(nargout(SS.function));
   end
%elseif isa(rFunc,'inline')
else
   Prob.OPTTB.r=rFunc;
   if strcmpi(rFunc,'ls2_rJS')
      % Now Tomlab calls Tomlab through intermediate fmincon interface - weird
      % Must set proper NumDiff, in case derivatives not defined
      if isempty(Prob.FUNCSX.J)
         funArgOut = 1;
      else
         funArgOut = 2;
      end
      funArgIn  = 2;
   else
      funArgOut = abs(nargout(rFunc));
      funArgIn  = abs(nargin(rFunc));
   end
   if funArgOut == 1
      Prob.NumDiff = max(1,Prob.NumDiff);
   else
      Prob.NumDiff=0;
   end
end

Prob.OPTTB.funArgOut=funArgOut;
Prob.OPTTB.funArgIn =funArgIn;

% Have to set Prob.LS.y for nlssolTL

if ~isfield(Prob.LS,'y')
   if funArgIn > 1
      r=feval(Prob.OPTTB.r,x_0,Prob.varargin{:});
   else
      r=feval(Prob.OPTTB.r,x_0);
   end
   Prob.LS.y = zeros(length(r),1);
else
   if isempty(Prob.LS.y)
      if funArgIn > 1
         r=feval(Prob.OPTTB.r,x_0,Prob.varargin{:});
      else
         r=feval(Prob.OPTTB.r,x_0);
      end
      Prob.LS.y = zeros(length(r),1);
   end
end
m = length(Prob.LS.y);

Prob.Solver.Alg  = DefPar(Prob.Solver,'Alg',0);
Solver           = DefPar(Prob.Solver,'Tomlab',[]);

if isempty(Solver)
   Solver = GetSolver(checkType([],Prob.probType),...
            (n > 500 | m > 2000) | Prob.LargeScale,0);
   Prob.Solver.Tomlab = Solver;
end

if Diagnostic
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
   fprintf(' Objective: TOMLAB ls_f and optim_rJ are calling  %s',rFunc);
   fprintf('\n');
   fprintf(' Gradient:  TOMLAB ls_g and optim_J are');
   if Prob.NumDiff==0
      fprintf(' using Jacobian in global variable LS_J\n');
   else
      fprintf(' using finite-differencing method # %d\n',Prob.NumDiff);
   end
   fprintf('\n');
   fprintf(' Number of lower bound constraints:              %d\n',...
           sum(~isinf(x_L)));
   fprintf(' Number of upper bound constraints:              %d\n',...
           sum(~isinf(x_U)));

   if isempty(Prob.b_U) | isempty(Prob.b_L)
      mEQ = 0;
   else
      mEQ=sum(~isinf(Prob.b_U) & Prob.b_L==Prob.b_U);
   end
   if isempty(Prob.b_U) 
      mU = 0;
   else
      mU =sum(~isinf(Prob.b_U))-mEQ;
   end
   if isempty(Prob.b_L) 
      mL = 0;
   else
      mL =sum(~isinf(Prob.b_L))-mEQ;
   end

   fprintf(' Number of linear equality constraints:          %d\n',mEQ);
   if mL > 0
      fprintf(' Number of linear lower bound constraints:       %d\n',mL);
   end
   fprintf(' Number of linear upper bound constraints:       %d\n',mU);
   fprintf('\n');
   fprintf('Algorithm selected is %s Alg # %d\n',Solver,Prob.Solver.Alg);
   fprintf('\n');
   disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
   fprintf(' End diagnostic information\n');
   fprintf('\n');
end

Result   = tomRun(Solver,Prob,Prob.PriLev);

x        = x_0;
x(:)     = Result.x_k;
f_k      = 2*Result.f_k;
r_k      = Result.r_k;
J_k      = Result.J_k;
ExitFlag = Result.ExitFlag;
Inform   = Result.Inform;

% Convert ExitFlag to OPTIM TB 
if ExitFlag==0
   ExitFlag=1;
   if PriLev > 0
      fprintf('lsqnonlin (%s): ',Solver);
      fprintf('Optimization Terminated successfully\n');
   end
else
   switch Inform
      case 101 
        ExitFlag=0;
        if PriLev > 0
           fprintf('lsqnonlin (%s): Too many iterations\n',Solver);
        end
      otherwise
        ExitFlag=-Inform;
        if PriLev > 0
           fprintf('lsqnonlin (%s): Possible Error.',Solver);
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

if nargout > 4
   Output.iterations  = Result.Iter;
   Output.funcCount   = n_r;
   Output.algorithm   = [Solver ': ' Result.SolverAlgorithm];
   if isempty(alphaV)
      Output.stepsize   = 1;
   else
      Output.stepsize   = alphaV(length(alphaV));
   end
end

if nargout > 5
   [v_k, Zv, P] = LagMult(Result.Prob,Result);

   ix = find(P(1:n) > 0);
   Lambda.upper      = zeros(n,1);
   Lambda.upper(ix)  = v_k(ix);

   ix = find(P(1:n) == -1);
   Lambda.lower      = zeros(n,1);
   Lambda.lower(ix)  = v_k(ix);
end

% MODIFICATION LOG
%
% 011204 hkh Test is nonempty before setting value from options into Prob.
% 011204 hkh Small changes in output.
% 030128 hkh Major revision for Tomlab v4.0
% 031201 hkh Change field AutoDiff to ADObj
% 040203 hkh Print ExitText if possible error (and PriLev > 0)
% 040728 med Pragma removed
% 041123 hkh Change call to tomRun
% 050422 hkh Handle new type of function handle in Matlab 7.x
% 050422 hkh Test and set NumDiff if calling TOMLAB solver
% 060814 med FUNCS used for callbacks instead
% 080414 hkh Revise help, add text about choice of solver
% 080521 med Defaults call added