% TOMLAB KNITRO (MI)NLP Solver
%
% function Result = knitroTL(Prob)
%
% INPUT:
%
% Prob   Problem structure in TOMLAB format
%
% --------------------------index.php---------------------
% Fields used in input structure Prob:
%
% x_0         Initial vector.
%
% x_L, x_U    Bounds on variables.
%
% A           Linear constraint matrix.
% b_L, b_U    Bounds on linear constraints.
%
% c_L, c_U    Bounds on nonlinear constraints.
%
% LargeScale  If 1, indicates that the problem should be treated as sparse.
%             This makes knitroTL send a sparse version of Prob.A to the
%             solver. To avoid poor performance, sparse ConsPattern and
%             d2LPattern should be given.
%             (Default 1)
%
% ConsPattern Sparsity pattern of the gradient of the nonlinear
%             constraints. Should be sparse if given.
%
% d2LPattern  Sparsity pattern of the Hessian of the Lagrangian function.
%
% f_Low       A lower bound on the objective function value. KNITRO will stop
%             if it finds an x for which f(x)<=f_Low.
%
% PriLevOpt   Print level in solver. TOMLAB /KNITRO supports solver
%             progress information being printed to the MATLAB command
%             window as well as to a file.
%
%             For output only to the command window: set PriLevOpt
%             to a positive value (1-6) and set:
%
%          >> Prob.KNITRO.PrintFile = '';
%
%             If PrintFile is set to a valid filename, the same
%             information will be written to the file (if PriLevOpt
%             is nonzero).
%
%             For printing only to the PrintFile, set a name (or
%             knitro.txt will be used) and a _negative_ PriLevOpt
%             value. For example:
%
%          >> Prob.KNITRO.PrintFile = 'myprob.txt';
%          >> Prob.PriLevOpt        = -2;
%
%             To run TOMLAB /KNITRO silently, set PriLevOpt=0; (default)
%
% MIP.IntVars
%             If empty, all variables are assumed non-integer
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer
%             variables
%
% -----------------------------------------------
% Fields used in Prob.KNITRO:
% -----------------------------------------------
%
% PrintFile   Name of file to print solver progress and status
%             information to. Please see PriLevOpt parameter
%             described above.
%
% objType     Flag telling the type of the objective function. TOMLAB can
%             determine this automatically for problems of type 'lp' and 'qp'
%             (see help checkType for more information on problem types)
%             but for a general user defined problem where the user has knowledge
%             of the objective function, this can be set to
%
%             0 - general/nonlinear objective
%             1 - linear objective function
%             2 - quadratic objective function
%
%             By passing this information to KNITRO, it can make certain
%             specializations for problems where all constraints are linear
%             and the objective is either linear or quadratic.
%
% mpec        List of complementary variables. If given directly by the user,
%             this must be a dense matrix with two columns containing the
%             indices of variables that should be complementary to each
%             other.
%
%             See mpecDemo and/or mcpAssign for a more flexible way to
%             specify complementarity properties. The quickguide may be
%             useful as well (lcpQG, qcpQG and mcpQG).
%
%
% checkderiv  Flag controlling the 1st derivative checking option before
%             solving. Allowed values:
%
%          1  Check first derivatives, then solve problem.
%          2  Only check derivatives, do not solve.
%
%             The tol_abs and tol_rel options control the absolute and
%             relative thresholds for what values to report. Derivative
%             elements with larger discrepancies than tol_abs and/or tol_rel
%             are printed.
%
% tol_abs     Absolute tolerance for derivative check. Default 1E-6.
%
% tol_rel     Absolute tolerance for derivative check. Default 1E-6.
%
%
% options     Structure with options values for KNITRO. The following
%             subfields in Prob.KNITRO.options may be set. See the manual
%             for a full list of parameters.
%
%  ALG        Choice of NLP algorithm
%
%             0 - Automatic selection (default)
%             1 - Barrier/Direct
%             2 - Barrier/CG
%             3 - Active set strategy
%
%  MIP_METHOD Choice of MINLP algorithm
%
%             0 - Automatic selection (default)
%             1 - Branch and bound
%             2 - Hybrid Quesada-Grossman
%
%  HESSOPT    Specifies how to calculate the Hessian of the Lagrangian function.
%
%             1 - KNITRO expects the user to provide the exact Hessian (default).
%
%             2 - KNITRO will compute a (dense) quasi-Newton BFGS Hessian.
%
%             3 - KNITRO will compute a (dense) quasi-Newton SR1 Hessian.
%
%             4 - KNITRO will compute Hessian-vector products using finite differences.
%
%             5 - KNITRO expects the user to provide a routine to
%                 compute the Hessian-vector products. If this option is
%                 selected, the calculation is handled by the Tomlab
%                 function "nlp_d2Lv.m".
%
%             6 - KNITRO will compute a limited-memory quasi-Newton BFGS.
%
%             NOTE: Typically if exact Hessians (or exact Hessian-
%             vector products) cannot be provided by the user but exact
%             gradients are provided and are not too expensive to
%             compute, option 4 above is recommended. The finite-
%             difference Hessian-vector option is comparable in terms
%             of robustness to the exact Hessian option (assuming exact
%             gradients are provided) and typically not too much slower
%             in terms of time if gradient evaluations are not the
%             dominant cost.
%
%             However, if exact gradients cannot be provided
%             (i.e. finite-differences are used for the first
%             derivatives), or gradient evaluations are expensive, it
%             is recommended to use one of the quasi-Newton options, in
%             the event that the exact Hessian is not
%             available. Options 2 and 3 are only recommended for small
%             problems (n < 1000) since they require working with a
%             dense Hessian approximation. Option 6 should be used in
%             the large-scale case.
%
%             NOTE: Options HESSOPT=4 and HESSOPT=5 are not available
%             when ALG=1.
%
% -----------------------------------------------------------------------
%
% Termination Message: At the end of the run a termination message is
% printed indicating whether or not the optimal solution was found and
% if not, or why the code terminated. See the manual for more information.
%
% -----------------------------------------------------------------------
%
% OUTPUT:
%
% Result   Structure with results (see ResultDef.m):
%
% x_k      Solution vector.
% x_0      Initial solution vector.
%
% f_k      Function value at optimum.
% g_k      Gradient of the objective function.
% H_k      Hessian of the Lagrangian function.
%
% c_k      Nonlinear constraint residuals.
% cJac     Nonlinear constraint gradients.
%
% xState   State of variables. Free == 0; On lower == 1; On upper == 2;
%          Fixed == 3;
%
% bState   State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
%
% cState   State of nonlinear constraints. Free == 0; Lower == 1; Upper == 2;
%          Equality == 3;
%
% v_k      Lagrangian multipliers (for bounds + dual solution vector).
%
% ExitFlag Exit status.
%
% Inform   KNITRO information parameter. See termination messages above.
%
% rc       Reduced costs. If ninf=0, last m == -v_k.
% Iter     Number of iterations.
% FuncEv   Number of function evaluations.
% GradEv   Number of gradient evaluations.
% ConstrEv Number of constraint evaluations.
% QP.B     Basis vector in TOMLAB QP standard.
% MinorIter Number of minor iterations.
% Solver   Name of the solver (knitro).
% SolverAlgorithm  Description of the solver.
%
% Result.KNITRO   Structure with KNITRO specifiec outputs. Fields set:
%
%   iv        Information vector:
%
%      iv(1)  Number of f(x) and c(x) evaluations
%      iv(2)  Number of g(x) and dc(x) evaluations
%      iv(3)  Number of Hessian evaluations
%      iv(4)  Number of Hessian*vector evaluations
%      iv(5)  Number of major iterations (also in Result.Iter)
%      iv(6)  Number of minor iterations (also in Result.MinorIter)
%      iv(7)  Absolute feasibility error
%      iv(8)  Absolute optimality error

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written May 5, 2003.  Last modified Aug 11, 2009.

function Result = knitroTL(Prob)

if nargin < 1, error('knitroTL needs the Prob structure as input'); end

Prob.solvType = 3; % NLP (CON) solver

KNITRO   = DefPar(Prob,'KNITRO',[]);
options = DefPar(KNITRO,'options',[]);

% Hessian information, if required by choice in options.HESSOPT
HL = triu(DefPar(Prob,'d2LPattern',[]));

if isfield(options,'gradopt')
   GRADOPT = options.gradopt;
elseif isfield(options,'GRADOPT')
   GRADOPT = options.GRADOPT;
else
   GRADOPT = 1;
end
if isempty(GRADOPT), GRADOPT = 1; end
if Prob.NumDiff == 6 & ~any(GRADOPT==[2 3])
   % Internal differentiation
   options.gradopt = 2;
   GRADOPT = 2;
end

if isfield(options,'hessopt')
   HESSOPT = options.hessopt;
elseif isfield(options,'HESSOPT')
   HESSOPT = options.HESSOPT;
else
   HESSOPT = [];
end
if isempty(HESSOPT)
   if (isempty(Prob.FUNCS.g) & Prob.ADObj <= 0) |  ...
         (isempty(Prob.FUNCS.dc) & Prob.mNonLin > 0 & Prob.ADCons <= 0)
      % The user has not provided any gradient or constraint Jacobian, and no
      % automatic differentiation is used to provide the 1st order information
      % Avoid two levels of differentiation, by default
      if Prob.LargeScale == 1
         if isempty(HL)
            % The problem is large, and the user has not provided the pattern
            % for the second derivative of the Lagrangian
            % Sparse BFGS will not need this information
            options.hessopt = 6;   % Sparse BFGS
            HESSOPT = 6;
         else
            % The problem is large, but the user has provided the pattern
            % for the second derivative of the Lagrangian
            % Therefor we are able to utilize sparse Hessian products
            options.hessopt = 4;   % Sparse Hessian vector products
            HESSOPT = 4;
         end
      else
         % The problem is NOT large. There is a choice between just using
         % 1st order information (Dense BFGS) or using finite differences
         % with sparse Hessian products. Seem reasonable to use Dense BFGS
         options.hessopt = 2;   % Dense BFGS
         HESSOPT = 2;
      end
   elseif (isempty(Prob.FUNCS.H) & Prob.ADObj >= 0) | ...
         (isempty(Prob.FUNCS.d2c) & Prob.mNonLin > 0 & Prob.ADCons >= 0)
      % The user has provided high accuracy gradient and constraint Jacobians
      % However, no explicit second order information or
      % Hessian/Lagrangians computed with automatic differentiation
      if Prob.LargeScale == 1
         if isempty(HL)
            % The problem is large, but the user has not provided the pattern
            % for the second derivative of the Lagrangian
            % Check the number of nonlinear constraints, if > 20
            if Prob.mNonLin > 20
               options.hessopt = 6;   % Sparse BFGS
               HESSOPT = 6;
            else
               % Low number of nonlinear constraints, try estimating
               % second order Sparse Hessian vector products
               options.hessopt = 4;   % Sparse Hessian vector products
               HESSOPT = 4;
            end
         else
            % The problem is large. The user has provided the pattern
            % for the second derivative of the Lagrangian
            options.hessopt = 4;   % Sparse Hessian vector products
            HESSOPT = 4;
         end
      else
         % Small problem, no patterns given
         % Use either a dense BFGS, avoiding use of second order information
         % Or use Hessian-vector finite differences (Ziena recommendation)
         % options.hessopt = 2;   % Dense BFGS
         options.hessopt = 4;     % Hessian-vector products, finite diffs
         HESSOPT = 4;
      end
   else
      % Two levels of derivative given, analytic or AD
      HESSOPT = 1;
   end
end

if abs(Prob.NumDiff) == 6 & ~any(HESSOPT==[2 3 4 6])
   % Internal differentiation?
   if Prob.LargeScale == 1
      if isempty(HL)
         options.hessopt = 6;   % Sparse BFGS
         HESSOPT = 6;
      else
         options.hessopt = 4;   % Sparse Hessian vector products
         HESSOPT = 4;
      end
   else
      options.hessopt = 2;   % Dense BFGS
      HESSOPT = 2;
   end
end

if any(HESSOPT == [1 5])
   % d2LPattern is used to decide on HESSOPT so have to set
   % it after this.
   if ~isempty(Prob.FUNCS.H) & strmatch(funch2str(Prob.FUNCS.H),{'lp_H','qp_H'},'exact')
      if isempty(Prob.FUNCS.c) & isempty(Prob.d2LPattern)
         if ~isempty(Prob.QP.F)
            Prob.d2LPattern = spones(Prob.QP.F);
         else
            Prob.d2LPattern = sparse([], [], [], Prob.N, Prob.N, 0);
         end
         Prob = iniSolve(Prob,3,0,0);
      else
         if strmatch(funch2str(Prob.FUNCS.H),'qp_H','exact')
            if isempty(Prob.HessPattern)
               Prob.HessPattern = spones(Prob.QP.F);
            end
         else
            if isempty(Prob.HessPattern)
               Prob.HessPattern = sparse([], [], [], Prob.N, Prob.N, 0);
            end
         end
         Prob = iniSolve(Prob,3,2,2);
      end
   else
      Prob = iniSolve(Prob,3,2,2);
   end
else
   Prob = iniSolve(Prob,3,1,1);
end

% Hessian information, if required by choice in options.HESSOPT
HL = triu(DefPar(Prob,'d2LPattern',[]));

Result=ResultDef(Prob);
Result.Solver='KNITRO';
PriLev=Prob.PriLevOpt;

if isempty(Prob.Name)
   Prob.Name = ['Problem ' num2str(Prob.P)];
end

BIG=1E20;

[bl, bu, n, m1, m2] = defblbu(Prob, BIG, 1);

m  = m1+m2;

% Variable bounds separately
xl = bl(1:n);
xu = bu(1:n);

cl = bl(n+1:end);
cu = bu(n+1:end);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n,  error('Linear constraints A MUST have n columns!'); end
   if mA~=m1, error('Linear constraints A MUST have m1 rows!'); end
end

% Logical vector for integers
IV = zeros(n,1);
IntVars=DefPar(Prob.MIP,'IntVars',[]);
if isempty(IntVars)
   mip = 0;
   IV = [];
elseif any(IntVars==0) | all(IntVars==1)
   % Assume binary logical vector given
   mip = 1;
   IV(1:length(IntVars)) = logical(IntVars);
else
   if any(IntVars < 1 | IntVars > n)
      error('knitro: Illegal IntVars vector');
   end
   IV(IntVars)=1;
   % Binary variables
   glcolidx = find(IV);
   ix = find(xl(glcolidx)==0 & (xu(glcolidx)==1) & (IV(glcolidx)==1));
   if ~isempty(ix)
      IV(glcolidx(ix)) = 2;
   end
   mip = 1;
end

% Safeguarded starting point
x_0 = DefPar(Prob,'x_0',zeros(n,1));
x_0 = max( bl(1:n),min(bu(1:n),x_0(:) ) );

Result.f_0 = nlp_f(x_0,Prob);
Result.x_0 = x_0;

% Lower bound on f(x)
f_Low = DefPar(Prob,'f_Low',-1E300);

LargeScale = DefPar(Prob,'LargeScale',1);

if LargeScale | issparse(Prob.A)
   Prob.A = sparse(Prob.A);
   nz = nnz(Prob.A);
else
   nz = mA*n;
end

if isempty(Prob.ConsPattern)
   Prob.ConsPattern = sparse( ones(m2,n) );
else
   Prob.ConsPattern = sparse( Prob.ConsPattern );
end
nz = nz + nnz(Prob.ConsPattern);

if m2 > 0
   % Determine the sparse problem structure
   if ~isempty(Prob.ConsPattern)
      [ix,iy]=find(Prob.ConsPattern);

      % Send linear index from multiple subscripts for nonzero pattern
      Prob.ConsIdx = sub2ind(size(Prob.ConsPattern),ix,iy);
   end
end

if isfield(options,'ALG')
   alg = options.ALG;
elseif isfield(options,'alg')
   alg = options.alg;
else
   alg = 0;
end
if isempty(alg)
   alg = 0;
end
if mip
   if isfield(options,'MIP_METHOD')
       alg = options.MIP_METHOD+4;
   else
       alg = 4;
   end
end

switch(alg)
   case 1,
      S='Interior/Direct NLP KNITRO';
   case 2,
      S='Interior/CG NLP KNITRO';
   case 3,
      S='Active Set SLQP NLP KNITRO';
   case 4,
      S='Default MINLP KNITRO';
   case 5,
      S='Branch and bound MINLP KNITRO';
   case 6,
      S='Hybrid Quesada-Grossman MINLP KNITRO';
   otherwise
      S='Default NLP KNITRO';
end
Result.SolverAlgorithm=S;

optParam = Prob.optParam;

% Max iters
if ~isfield(options,'MAXIT') & ~isfield(options,'maxit')
   options.MAXIT = DefPar(optParam,'MaxIter',2000);
end

% Max time
if ~isfield(options,'MAXTIMECPU') & ~isfield(options,'maxtimecpu')
   options.MAXTIMECPU = DefPar(Prob,'MaxCPU',1e8);
end

% Type of objective function
objType = DefPar(KNITRO,'objType',[]);

% Detect LP or QP from Prob.probType but only if the user has not already set this.
if isempty(objType)
   if checkType('qp',Prob.probType)==1
      objType = 2;
   elseif checkType('lp',Prob.probType)==1
      objType = 1;
   else
      objType = 0;
   end
elseif ~any( [0,1,2]==objType )
   error('Illegal Prob.KNITRO.objType setting - 0,1,2 allowed');
end

% Lin/quad objective only taken into account if all constraints are linear
if m2>0
   objType = 0;
end

% KNITRO callback to m-file?
cb=DefPar(KNITRO,'callback','');
if ~isempty(cb)
   % SAL users - change this!
   switch(exist(cb,'file'))
      case {2,3,6}
         % m, mex, or p-file's are ok. Could check # of in/out too.
      otherwise
         warning(['The callback function given (''' cb ''') is not an m-, mex-, or p-file. Callbacks will be disabled']);
         cb = '';
   end
end

% PrintFile
PrintFile = DefPar(KNITRO,'PrintFile',[]);

% Negative printlevel means we want printing to a file only. If no
% name is given, set knitro.txt as default name.
if PriLev < 0 & isempty(PrintFile)
   PrintFile = 'knitro.txt';
end

% No Hessian d2L(x) information available? We need this only for HESSOPT=1
if HESSOPT==1
   if isempty(HL)

      % For LargeScale problems we get this by evaluating the second
      % derivatives at the starting point. Not 100% safe though.

      if LargeScale
         if ~isempty(Prob.HessPattern)
            HL = spones(Prob.HessPattern);
         else
            HL = estHessPattern(Prob);
         end
         if m2>0
            if isfield(Prob,'d2cPattern')
               if ~isempty(Prob.d2cPattern)
                  d2c = spones(Prob.d2cPattern);
               else
                  d2c = estd2cPattern(Prob);
                  Prob.d2cPattern = d2c;
               end
            else
               d2c = estd2cPattern(Prob);
               Prob.d2cPattern = d2c;
            end
            HL = spones(HL + d2c);
         end
      else
         % Non-large scale - dense Hessian.
         HL = sparse(triu(ones(n,n)));
      end
   end

   HL = triu(HL);

   if ~issparse(HL), HL = sparse(HL); end
else
   HL = [];
end

% mpec pairs
mpec = DefPar(KNITRO,'mpec',[]);
if ~isempty(mpec)
   if size(mpec,2) ~= 2
      error('knitroTL: Prob.KNITRO.mpec must have exactly two columns')
   end
   if sum(any(mpec>Prob.N)) | sum(any(mpec<1))
      error('knitroTL: Prob.KNITRO.mpec contains illegal values < 1 or > Prob.N')
   end
   if issparse(mpec)
     mpec = full(mpec);
   end
end

% TODO: implement ctype. The mex does a rudimentary check; linear
% constraints = convex and nonlinear equalities = nonconvex
ctype = [];

try
% Call the solver
[Inform,iv,x_k,f_k,g_k,c_k,v_k,cJac,H_k] = ...
   Tknitro(n,m1,m2,...
   xl,xu,x_0,...
   Prob.A,cl,cu,...
   Prob.ConsPattern,nz,HL,...
   mpec, objType, f_Low,...
   PriLev,PrintFile,options,...
   cb,Prob,IV,ctype);
catch
   l=lasterror;
   if(strcmp(l.identifier,'MATLAB:invalidMEXFile'))
      tomlabsharederror;
   else
      rethrow(l);
   end
end

% Tknitro version 4.0 returns H_k only for HESSOPT = 1;
if ~isempty(H_k)
   H_k = triu(H_k,1) + diag(diag(H_k)) + triu(H_k,1)';
end

Result.Iter      = iv(5);
Result.MinorIter = iv(6);
Result.x_k  = x_k;
Result.x_0  = x_0;
Result.f_k  = f_k;
Result.g_k  = g_k;
Result.H_k  = H_k;
% KNITRO uses negative lagrange multipliers. Negate to
% convert to TOMLAB standard.
Result.v_k  = -[ v_k(m+1:end) ; v_k(1:m) ];
Result.cJac = cJac;

% Separate linear and nonlinear constraints
% if m1>0, Result.Ax  = c_k(1:m1);   else Result.Ax = [];  end
% if m2>0, Result.c_k = c_k(m1+1:m); else Result.c_k = []; end

if m1>0, Result.Ax  = Prob.A*x_k; else Result.Ax = []; end
if m2>0, Result.c_k = nlp_c(x_k,Prob); else Result.c_k = []; end

% StateDef wants bl,bu to be [vars,nonlin,lin], unlike what is returned by
% defblbu above, therefore set additional input Order to 1 to reverse order.

% The barrier method will put x variables slightly outside bounds
FEASTOL = DefPar(options,'FEASTOL',1E-6);

Result = StateDef(Result, x_k(1:n), Result.Ax, Result.c_k, ...
   FEASTOL, optParam.bTol, optParam.cTol, bl, bu, 1);

% Inform -> ExitFlag, ExitText translation
switch(Inform)
   % Success
   case 0
      ExitText = 'Locally optimal solution found';
   case -100
      ExitText = 'Primal feasible solution estimate cannot be improved. It appears to be optimal, but desired accuracy in dual feasibility could not be achieved.';
   case -101
      ExitText = 'Primal feasible solution; terminate because the relative change in solution estimate < xtol.';
   case -102
      ExitText = 'Primal feasible solution estimate cannot be improved; desired accuracy in dual feasibility could not be achieved.';

      % Infeasibility
   case -200
      ExitText = 'Convergence to an infeasible point. Problem may be locally infeasible.';
   case -201
      ExitText = 'Terminate at infeasible point because the relative change in solution estimate < xtol.';
   case -202
      ExitText = 'Current infeasible solution estimate cannot be improved. Problem may be badly scaled or perhaps infeasible.';
   case -203
      ExitText = 'MULTISTART: No primal feasible point found.';

      % Unboundedness
   case -300
      ExitText = 'Problem appears to be unbounded. Iterate is feasible and objective magnitude > objrange.';

      % Resource limit
   case -400
      ExitText = 'Iteration limit reached.';
   case -401
      ExitText = 'Time limit reached.';
   case -403
      ExitText = 'All nodes have been explored.';
   case -404
      ExitText = 'Terminating at first integer feasible point.';
   case -405
      ExitText = 'Subproblem solve limit reached.';
   case -406
      ExitText = 'Node limit reached.';
   case -499
      ExitText = 'Terminated by the user.';

    % Other error conditions
   case -500
      ExitText = 'Callback function error.';
   case -501
      ExitText = 'LP solver error.';
   case -502
      ExitText = 'Evaluation error.';
   case -503
      ExitText = 'Not enough memory available to solve problem.';

   otherwise
   if Inform >= -505 && Inform <= 600
      ExitText = 'Input error. Run with PriLevOpt>0 for details';
   else
      ExitText = 'Unknown error';
   end
end

if Inform > -200
   ExitFlag = 0;
elseif Inform > -300 && Inform <= -200
   ExitFlag = 4; % Infeasible
elseif Inform > -400 && Inform <= -300
   ExitFlag = 2; % Unbounded
elseif Inform > -500 && Inform <= -400
   ExitFlag = 1; % Time or iter limit
else
   ExitFlag = 10; % Input error in some form
end


Result.ExitText = ExitText;
Result.ExitFlag = ExitFlag;
Result.Inform   = Inform;

Result=endSolve(Prob,Result);

% MODIFICATION LOG
%
% 030505 ango Wrote file
% 030708 ango Hessian handling added
% 030729 ango Comments update
% 030829 ango Further comments updates
% 030903 ango ExitFlag/Inform handling
% 030917 ango f_Low implemented; LargeScale=1 is default.
% 031211 ango Changed the name of this file from TknitroTL
% 040102 hkh  Revision for v4.2, call iniSolve and endSolve
% 040109 hkh  Skip call to ProbCheck, use returned f_k and g_k
% 040203 ango LargeScale Hessian handling modified.
% 040602 med  Help fix PriLev to PriLevOpt
% 041005 ango ISLP, ISQP checked and set for LP, QP problems.
% 041110 ango Help updated for KNITRO 4, MaxCPU added
% 041111 med  Help fixes, MAXIT and MAXTIME fixed.
% 041128 hkh  Different calls to iniSolve dependent on HESSOPT
% 041129 med  optPar settings removed
% 041202 hkh  Removed inefficient call to Prob.optParam, use optParam instead
% 041202 hkh  Use new StateDef (and defblbu), avoid reversing bl/bu order
% 041202 hkh  Use FEASTOL instead of xTol, interior point methods are inexact
% 041210 hkh  Revise GRADOPT,HESSOPT use, avoid two level of differentiation
% 050223 frhe Lagrange multipliers negated
% 050421 hkh  Revise default HESSOPT use, check 2nd order information
% 050421 hkh  Separate treatment of A and ConsPattern, avoid copying
% 050804 ango Changes to Hessian pattern generation for HESSOPT=1
% 050902 ango Add a triu(HL), slight hessian pattern change
% 050902 med  Moved HL to after iniSolve call (and keep before)
% 050908 med  d2LPattern set automatically for LP/QP/LPCON/QPCON problems before iniSolve call
% 051216 med  Help updated
% 060524 med  MPEC added to help
% 060607 ango Change MPEC help and add checks on Prob.KNITRO.mpec size & values
% 060818 med  FUNCS used for callbacks
% 061004 med  Help updated (MAXTIME)
% 061026 ango Remove obsolete GRADOPT = 4 and 5
% 061212 med  MPEC help updated
% 061213 med  Shared error added
% 070307 med  Updated for version 5.1.1
% 070823 ango HESSOPT=1 revised, no need to get row/col data
% 070830 ango Always get ConsPattern - the MEX needs this
% 070905 med  Fixed lp_H and qp_H checks for function handle
% 080507 ango Update for 5.2.0 (Inform values changed)
% 080929 ango Fix error in callback feature
% 081030 ango Gradient check
% 090226 ango KNITRO 6 updates
% 090804 med  Removed version check
% 090811 med  Termination messages removed, see manual
