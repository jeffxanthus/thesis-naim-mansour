% slsSolve.m:
%
% Finds a Sparse Least Squares (sls) solution to a constrained least squares
% problem with many variables, with the use of a suitable NLP TOMLAB solver
%
% The use is mainly for sparse problems, where the structure in both the
% Jacobian of the residuals and the Jacobian of the nonlinear constraints
% are utilized by a sparse NLP solver, e.g. SNOPT.
%
% function Result = slsSolve(Prob, PriLev)
%
% Minimization problem:
%
%        min  0.5*r'r, where r(x) is in R^m
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% The constrained least squares (L2) problem is solved in slsSolve by
% rewriting the problem as a general constrained optimization problem.
% A set of m additional variables z, stored as x(n+1:m), is added
%
%        min    0.5*z'z,
%         x
%        s/t   x_L <=   x              <= x_U
%              b_L <= A x              <= b_U
%              c_L <= c(x)             <= c_U
%                0 <= r(x) - z         <= 0
%
% ---------------------------------------------------------------------------
%
%
% INPUT PARAMETERS
%   Prob       Structure Prob. Prob must be defined.
%              Best is to use Prob = clsAssign(.....), if using the TQ format.
%              The problem should be created in the TOMLAB constrained
%              nonlinear least squares format (cls)
%              This format is aimed for a vector valued function, with a
%              Jacobian matrix as the derivative.
%
%   PriLev     The second input argument. Default == 2.
%              If PriLev == 0, slsSolve is silent, except for error messages.
%              If > 0, slsSolve prints summary output about problem
%              transformation
%              slsSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%              PriLev == 2 displays standard output in PrintResult.
%
%   Use Prob = probInit('name of file',problem_number'); if solving
%   a predefined problem in the Init File (IF) format.
%
%   Extra fields used in Prob:
%
%   SolverL2    Name of the TOMLAB solver. Valid names are:
%               conSolve, nlpSolve, (sTrustr, clsSolve)
%               If TOMLAB /SOL is installed: minos, snopt, npsol.
%               optPar(66) is set to 1 when running snopt7.
%
%   L2Type     =1 standard constrained formulation
%
%
%   The rest of the fields in Prob should be defined as wanted by the
%   selected solver. See the help for the solver.
%
%   In particular:
%
%   x_0:     Starting point
%   x_L:     Lower bounds for x
%   x_U:     Upper bounds for x
%   b_L:     Lower bounds for linear constraints
%   b_U:     Upper bounds for linear constraints
%   A:       Linear constraint matrix
%   c_L:     Lower bounds for nonlinear constraints
%   c_U:     Upper bounds for nonlinear constraints
%
%   ConsPattern The pattern of the constraint Jacobian
%   JacPattern  The pattern of the residual Jacobian
%   Note that Prob.LS.y must have the correct residual length, if JacPattern
%   is empty (but not ConsPattern).
%   slsSolve will create the new Prob.ConsPattern to be used by the solver
%   using the information in ConsPattern and JacPattern.
%
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.r_k, Result.J_k,
%        Result.c_k, Result.cJac, Result.x_0, Result.xState, Result.cState,
%        Result.v_k, is transformed back to the original problem.
%        Result.g_k is Result.J_k'*Result.r_k.
%        The output in Result.Prob is the result after slsSolve
%        transformed the problem, i.e. the altered Prob structure
%
% EXAMPLE:
%
% See slsDemo.m in tomlab\examples

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Apr 13, 2002.   Last modified July 24, 2011.

function Result = slsSolve(Prob, PriLev)

%#function L2_f L2_g L2_H L2_c L2_dc

if nargin < 2
   PriLev = [];
   if nargin < 1
      error('slsSolve needs input structure Prob');
   end
end
if isempty(PriLev), PriLev = 2; end

if isfield(Prob,'SolverL2')
   Solver=Prob.SolverL2;
   Solver=deblank(Solver);
else
   Solver = [];
end
if isfield(Prob,'L2Type')
   L2Type=Prob.L2Type;
else
   L2Type = [];
end
if isempty(L2Type), L2Type=1; end

solvType=checkType('con');

if ~isempty(Prob.JacPattern),  Prob.LargeScale = 1; end
if ~isempty(Prob.ConsPattern), Prob.LargeScale = 1; end

if isempty(Solver), Solver=GetSolver('con',1,0); end

Prob.CHECK = 0; % Force a test in ProbCheck
Prob=ProbCheck(Prob,Solver,solvType);

Prob.L2.r  = Prob.FUNCS.r;
Prob.L2.J  = Prob.FUNCS.J;
Prob.L2.c  = Prob.FUNCS.c;
Prob.L2.dc = Prob.FUNCS.dc;
args = zeros(9,1);
if ~isempty(Prob.FUNCS.r)
   args(7)    = xnargin(Prob.FUNCS.r);
end
if ~isempty(Prob.FUNCS.J)
   args(8)    = xnargin(Prob.FUNCS.J);
end
if ~isempty(Prob.FUNCS.c)
   args(4)    = xnargin(Prob.FUNCS.c);
end
if ~isempty(Prob.FUNCS.dc)
   args(5)    = xnargin(Prob.FUNCS.dc);
end

if L2Type >= 1
   % Always use now
   % Standard approach - L2Type = 1
   
   n = Prob.N;
   Prob = iniSolve(Prob,3,1,1);

   % Safe-guard starting point
   Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

   % Check the residual and Jacobian supplied
   Prob.Mode   = 2;
   Prob.nState = 1;
   r  = nlp_r(Prob.x_0,Prob);
   m  = length(r);
   Prob.L2.m    = m;
   if m == 0
      fprintf('Number of residuals supplied is 0\n');
      error('Illegal size of initial residual vector');
   end
   mmY = length(Prob.LS.y);
   if mmY == 0
      Prob.LS.y = zeros(m,1);
   elseif mmY ~= m
      fprintf('Number of residual elements is %d\n',m);
      fprintf('The number of observations in Prob.LS.y is %d\n',mmY);
      error('Should be equal, illegal size of vector Prob.LS.y');
   else
      Prob.LS.y = Prob.LS.y(:);
   end
   Prob.Mode   = 1;
   NumDiff     = Prob.NumDiff;
    ADObj       = Prob.ADObj;
   if NumDiff ~= 6
      J  = nlp_J(Prob.x_0,Prob);
      [m1,m2] = size(J);
      if m1 ~= m
         fprintf('Number of residual elements is %d\n',m);
         fprintf('The Jacobian of the residuals has %d rows\n',m1);
         error('Should be equal, illegal size of Jacobian');
      end
      if m2 ~= n
         fprintf('Number of variables is %d\n',n);
         fprintf('The constraint jacobian has %d columns\n',m2);
         error('Should be equal, illegal size of Jacobian');
      end
   end

   if ~isempty(Prob.A)
      mA = size(Prob.A,1);
      % Adjust A matrix for reformulated problem
      Prob.A = [sparse(Prob.A),sparse(mA,m)];
   else
      mA = 0;
   end

   % Check and adjust bounds on nonlinear constraints
   ML       = length(Prob.c_L);
   MU       = length(Prob.c_U);
   M        = max(ML,MU);
   [M1,M2]  =  size(Prob.ConsPattern);
   ConsDiff = Prob.ConsDiff;
   ADCons   = Prob.ADCons;
   if M > 0 | M1 > 0
      Prob.Mode = 2;
      c  = nlp_c(Prob.x_0,Prob);
      M  = length(c);
      if ML ~= M
         fprintf('Number of nonlinear constraints is %d\n',M);
         fprintf('Prob.c_L has %d rows\n',ML);
         error('Should be equal, illegal size of Prob.c_L');
      end
      if MU ~= M
         fprintf('Number of nonlinear constraints is %d\n',M);
         fprintf('Prob.c_U has %d rows\n',MU);
         error('Should be equal, illegal size of Prob.c_U');
      end
      Prob.Mode = 1;
      if ConsDiff ~= 6
         dc = nlp_dc(Prob.x_0,Prob);
         if ~isempty(dc)
            [Mc1,Mc2] = size(dc);
            if Mc1 ~= M
               fprintf('Number of nonlinear constraints is %d\n',M);
               fprintf('The constraint jacobian has %d rows\n',Mc1);
               error('Should be equal, illegal size of constraint Jacobian');
            end
            if Mc2 ~= n
               fprintf('Number of variables is %d\n',n);
               fprintf('The constraint jacobian has %d columns\n',Mc2);
               error('Should be equal, illegal size of constraint Jacobian');
            end
         end
      end
   end
   if ~isempty(Prob.JacPattern)
      [mJ1,mJ2] = size(Prob.JacPattern);
      if mJ1 ~= m
         fprintf('Number of residual elements is %d\n',m);
         fprintf('There are %d rows in Prob.JacPattern\n',mJ1);
         error('Should be equal, illegal size of Prob.JacPattern');
      end
      if mJ2 ~= n
         fprintf('Number of variables is %d\n',n);
         fprintf('Prob.JacPattern has %d columns\n',mJ2);
         error('Should be equal, illegal size of Prob.JacPattern');
      end
   end
   
   if ~isempty(Prob.ConsPattern)
      if M1 ~= M
         fprintf('Length of nonlinear constraint bounds is %d\n',M);
         fprintf('There are %d rows in Prob.ConsPattern\n',M1);
         error('Should be equal, illegal size of Prob.ConsPattern');
      end
      if M2 ~= n
         fprintf('Number of variables is %d\n',n);
         fprintf('Prob.ConsPattern has %d columns\n',M2);
         error('Should be equal, illegal size of Prob.ConsPattern');
      end
      if ~isempty(Prob.JacPattern)
         %m = size(Prob.JacPattern,1);
         Prob.ConsPattern = [sparse(Prob.ConsPattern),sparse(M,m); ...
              sparse(Prob.JacPattern),speye(m,m)];
      else
         m = length(Prob.LS.y);
         Prob.ConsPattern = [sparse(Prob.ConsPattern),sparse(M,m); ...
              sparse(ones(m,n)),speye(m,m)];
      end
   elseif ~isempty(Prob.JacPattern)
      %m = size(Prob.JacPattern,1);
      if M > 0
         Prob.ConsPattern = [sparse(ones(M,n)),sparse(M,m); ...
             sparse(Prob.JacPattern),speye(m,m)];
      else
         Prob.ConsPattern = [sparse(Prob.JacPattern),speye(m,m)];
      end
   end

   % Reformulate problem
   Prob.JacPattern = [];

   % Add m extra variables
   Prob.x_0 = [Prob.x_0(:);r];
   Prob.x_L = [Prob.x_L(:);-Inf*ones(m,1)];
   Prob.x_U = [Prob.x_U(:); Inf*ones(m,1)];
   Prob.x_min = [];
   Prob.x_max = [];

   Prob.N       = n + m;
   
   % Add m extra constraints from the m residuals
   Prob.c_L = [Prob.c_L(:);zeros(m,1)];
   Prob.c_U = [Prob.c_U(:);zeros(m,1)];
   Prob.mNonLin = length(Prob.c_L);

   if PriLev > 0
      fprintf('The problem has %d variables (Columns),',n);
      fprintf(' slsSolve adds unbounded variables (Columns) %d to %d\n', ...
                n+1,n+m);
      fprintf('These extra variables is the objective residual values\n');
      fprintf('\n');
      if mA == 0
         fprintf('The problem has %d linear constraints (Rows)\n',mA);
      else
         fprintf('The problem has %d linear constraints',mA);
         fprintf(' defined as constraint (Row) %d to %d\n',1,mA);
      end
      if M == 0
         fprintf('The problem has %d nonlinear constraints (Rows)\n',M);
      else
         fprintf('The problem has %d nonlinear constraints',M);
         fprintf(' defined as constraint (Row) %d to %d\n',mA+1,mA+M);
      end
      fprintf('The problem has %d residuals, by slsSolve ',m);
      fprintf('defined as constraint (Row) %d to %d',mA+M+1,mA+M+m);
      fprintf('\n');
   end

   Prob = tomFiles(Prob,'L2_f','L2_g','L2_H','L2_c','L2_dc');

   Prob.nState      = 0;
   Prob.L2.args     = args;
   Prob.L2.NumDiff  = NumDiff;
   Prob.L2.ConsDiff = ConsDiff;
   Prob.L2.ADObj    = ADObj;
   Prob.L2.ADCons   = ADCons;
   Prob.NumDiff     = 0;
   Prob.ConsDiff    = 0;
   Prob.ADObj       = 0;
   Prob.ADCons      = 0;
   if ~isempty(findstr('snopt', lower(Solver))) & Prob.SOL.optPar(66) == -999
       Prob.SOL.optPar(66) = 1;
   end
   Result=tomRun(Solver,Prob,PriLev);

   Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation
   % Return only original number of variables, remove f(x) value slack 
   Result.x_0     = Result.x_0(1:n);
   Result.xState  = Result.xState(1:n);
   Result.cState  = Result.cState(1:M);
   if ~isempty(Result.v_k)
      Result.v_k  = Result.v_k([1:n,n+m+1:n+m+mA+M]);
   end
   if ~isempty(Result.x_k)
      Result.r_k  = Result.x_k(n+1:n+m);
      Result.x_k  = Result.x_k(1:n);
   end
   if ~isempty(Result.c_k)
      Result.c_k  = Result.c_k(1:M);
   end
   if ~isempty(Result.cJac)
      Result.J_k  = Result.cJac(M+1:M+m,1:n);
      Result.cJac = Result.cJac(1:M,1:n);
      Result.g_k  = Result.J_k'*Result.r_k;
   else
      Result.g_k  = [];
   end
end

% MODIFICATION LOG:
%
% 020413  hkh  Written
% 020416  hkh  Fix handling of numerical derivatives, send Cons-/NumDiff
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin
% 040728  med  tomFiles used instead
% 041123  hkh  Change call to tomRun
% 041130  hkh  Prob.ConsDiff must be 0, not max(NumDiff,ConsDiff);
% 041201  hkh  Force a check in ProbCheck, setting Prob.CHECK=0
% 041222  med  Safeguard added for x_0
% 050421  hkh  Comments on some unnecessary tests
% 050422  hkh  Avoid Lagrange multipliers in PrintResult
% 050726  med  optPar(66) = 1 by default when using snopt7
% 050901  med  Removed size checks on x_L, x_U, x_0, A
% 050901  med  All pattern input now sparse
% 050901  med  Removed unnecessary variables
% 050901  med  Removed setting mLin, number of constr not changed
% 060814  med  FUNCS used for callbacks instead
% 070202  med  GetSolver LargeScale set to 1 at all times
% 110722  hkh  Numeric & AD differentation for both residual and constraints
% 110724  hkh  Test ConsDiff,NumDiff for ~=6, not < 6

