% infSolve.m:
%
% Finds a constrained weighted minimax solution of a function of several variables 
% with the use of any suitable TOMLAB nonlinear programming solver
%
% function Result = infSolve(Prob, PriLev)
%
% Minimization problem (without residual weighting):
%
%        min  max r(x), where r(x) is in R^m
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% The minimax problem is solved in infSolve by rewriting the problem as
% a general constrained optimization problem.
% One additional variable z, stored as x(n+1), is added
%
%        min    z,    where r(x) is in R^m
%         x
%        s/t   x_L <=   x(1:n)         <= x_U
%             -Inf <=   z              <= Inf
%              b_L <= A x              <= b_U
%              c_L <= c(x)             <= c_U
%             -Inf <= r(x) - z*e       <= 0,  e is in R^m, e(i)=1 for all i
%
% To handle cases where an element i in r(x), r_i(x) is taken the absolute
% value of:  min max |r_i(x)|, expand the problem with extra residuals
% with the opposite sign: [r_i(x); -r_i(x)]
%
% Residual weighting may be applied in a general way, see help clsAssign
% and the parameters weightType and weightY.
%
% ---------------------------------------------------------------------------
%
% INPUT PARAMETERS
%   Prob       Structure Prob. Prob must be defined.
%              Best is to use Prob = clsAssign(.....), if using the TQ format.
%              The problem should be created in the TOMLAB constrained
%              nonlinear least squares format (cls)
%              This format is aimed for a vector valued function, with a
%              Jacobian matrix as the derivative.
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%           Defines mixed-integer optimization parameters. Fields used:
%  IntVars:  
%           If empty, all variables are assumed non-integer 
%           If islogical(IntVars) (=all elements are 0/1), then
%           1 = integer variable, 0 = continuous variable.
%           If any element >1, IntVars is the indices for integer variables
%
%
%   PriLev     The second input argument. Default == 2.
%              If PriLev == 0, infSolve is silent, except for error messages.
%              If > 0, infSolve prints summary output about problem
%              transformation
%              infSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%              PriLev == 2 displays standard output in PrintResult.
%
%   Use Prob = probInit('name of file',problem_number'); if solving
%   a predefined problem in the Init File (IF) format.
%
%   Extra fields used in Prob:
%
%   SolverInf   Name of the TOMLAB solver. Valid names are:
%               conSolve (clsSolve)
%               If TOMLAB /SOL or /SNOPT is installed: minos, snopt
%               If TOMLAB /SOL or /NPSOL is installed: npsol (nlssol)
%               If TOMLAB /KNITRO is installed: KNITRO
%               If TOMLAB /CONOPT is installed: CONOPT
%               If clsSolve or NLSSOL is used, a LS approach is used
%
%   f_Low       A lower bound on the optimal function value.
%               Not crucial, if not set default -1E300 is used
%   f_Upp       An upper bound on the optimal function value.
%               Not crucial, if not set default 1E300 is used
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
%   ConsPattern The pattern of the constraint Jacobian (derivatives of c(x))
%   JacPattern  The pattern of the residual Jacobian (derivatives of r(x))
%
%   infSolve will create the new Prob.ConsPattern to be used by the solver
%   using the information in ConsPattern and JacPattern.
%
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.r_k, Result.J_k,
%        Result.c_k, Result.cJac, Result.x_0, Result.xState, Result.cState,
%        Result.v_k, is transformed back to the original problem.
%        Result.g_k is Result.J_k'*Result.r_k.
%        The output in Result.Prob is the result after infSolve
%        transformed the problem, i.e. the altered Prob structure
%
% EXAMPLE:
%
% See minimaxDemo.m in tomlab\examples

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written May 22, 1999.   Last modified July 24, 2011.

%   InfType     =1 constrained formulation. = 2 LS penalty approach (NOT OK)
% ---------------------------------------------------------------------------
% NOTE! The following approach is experimental. Use InfType == 1 !!!
% Another alternative problem formulation in infSolve is to solve the problem
% as a nonlinear least squares problem:
% One additional variable z, stored as x(n+1), is added
% Also m additional slack variables s, stored as x(n+2:n+m+1)
%
%        min      z^2 , where r(x) is in R^m    (InfType == 2)
%         x
%        s/t   x_L <=   x(1:n)       <= x_U, x in R^n
%             -Inf <=   z            <= Inf, z in R
%                0 <=   s            <= Inf, s in R^m, s slack variables
%              b_L <= A x            <= b_U
%              c_L <= c(x)           <= c_U
%             -Inf <= r(x) - z*e + s  = 0,  e is in R^m, e(i)=1 for all i
%
% The problem is solved with a penalty approach
% (tau increasing penalty parameter):
%
%        min      z'z + tau^2*(r-z*e + s)'*(r-z*e + s)
%         x
%        s/t   x_L <=   x(1:n)       <= x_U, x in R^n
%             -Inf <=   z            <= Inf, z in R
%                0 <=   s            <= Inf, s in R^m, s slack variables
%              b_L <= A x            <= b_U
%              c_L <= c(x)           <= c_U
%
% The full residual is formed as [tau*(r-z*e+s);z]
%
% The Jacobian structure is
%
%  J(x,z,s) = [J(x) -1*e I ]
%
% which is utilized

function Result = infSolve(Prob, PriLev)

%#function mmx_f mmx_g mmx_H mmx_c mmx_dc

if nargin < 2
   PriLev = [];
   if nargin < 1
      error('infSolve needs input structure Prob');
   end
end
if isempty(PriLev), PriLev = 2; end

if isfield(Prob,'SolverInf')
   Solver=Prob.SolverInf;
   Solver=deblank(Solver);
else
   Solver = [];
end
%if isfield(Prob,'InfType')
%   InfType=Prob.InfType;
%else
%   InfType = [];
%end
%if isempty(InfType), InfType=1; end

%if InfType == 2
%   solvType=checkType('cls');
%   if isempty(Solver), Solver=GetSolver('cls',Prob.LargeScale,0); end
%else
   solvType=checkType('con');
   if ~isempty(Prob.JacPattern),  Prob.LargeScale = 1; end
   if ~isempty(Prob.ConsPattern), Prob.LargeScale = 1; end
   if isempty(Solver), Solver=GetSolver('con',Prob.LargeScale,0); end
%end

Prob.CHECK = 0; % Force a test in ProbCheck
Prob=ProbCheck(Prob,Solver,solvType);

Prob.minimax.r  = Prob.FUNCS.r;
Prob.minimax.J  = Prob.FUNCS.J;
Prob.minimax.c  = Prob.FUNCS.c;
Prob.minimax.dc = Prob.FUNCS.dc;

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

%if InfType == 2 
%   if isempty(Prob.x_0)
%      n        = max([length(Prob.x_L),length(Prob.x_U),size(Prob.A,2)]);
%      Prob.x_0 = zeros(n,1);
%   else
%      n        = length(Prob.x_0);
%   end
%   r    = feval(Prob.FUNCS.r,Prob.x_0,Prob);
%   m    = length(r);
%   rMax = max(abs(r));
%
%   rho=10;
%
%   %z=rMax/2;
%   z=rMax;
%   s        = max(0,z-abs(r));
%   Prob.x_0 = [Prob.x_0;z;s];
%   Prob.N   = length(Prob.x_0);
%   Prob.x_L = [Prob.x_L;0;zeros(m,1)];
%   Prob.x_U = [Prob.x_U;rMax;rMax*ones(m,1)];
%   if ~isempty(Prob.A)
%      Prob.A = sparse([Prob.A,zeros(size(Prob.A,1),m+1)]);
%   end
%   Prob.mLin = size(Prob.A,1);
%
%   Prob.r_k = [r-z*ones(m,1)+s;z];
%
%   Prob     = tomFiles(Prob,'ls_f','ls_g','ls_H', ...
%                     Prob.FUNCS.c,Prob.FUNCS.dc,Prob.FUNCS.d2c, 'mmx_r','mmx_J');
%
%   %Prob.optParam.MaxIter=100;
%
%   epsX                 = Prob.optParam.eps_x;
%   epsG                 = Prob.optParam.eps_g;
%   Prob.optParam.eps_x  = max(1E-3,epsX);
%   Prob.optParam.eps_g  = 1E-4;
%   Prob.LS.weightY    = [rho*ones(m,1);1];
%   Prob.LS.weightType = 1;
%
%   Result=tomRun(Solver,Prob,2);
%
%   fprintf('\n');
%   xprinte(Result.x_k,'x_k:');
%
%   %fprintf('ResEv %d - JacEv %d - Iter %d\n',n_r,n_J,Result.Iter);
%   fprintf('Minimax L2-result %30.20f\n',max(Result.r_k));
%
%   xOld =Result.x_k;
%   notConv = 1
%   Iter = 1;
%
%   while notConv & Iter < 12
%      rho = rho * 10;
%      %Prob.LS.weightY(end) = rho;
%      Prob.LS.weightY    = [rho*ones(m,1);1];
%
%      % Weighted NLLS step
%      %Prob.LS.weightY=abs(r).*Prob.LS.weightY;
%      %wSum=sum(Prob.LS.weightY);
%      %Prob.LS.weightY=Prob.LS.weightY/wSum;
%      %wY=Prob.LS.weightY;
%      %xprinte(wY,'wY:');
%      %%wSum=sum(Prob.LS.weightY)
%
%      Prob.x_0 = xOld;
%
%      Result=tomRun(Solver,Prob,2);
%
%      ResEv=ResEv + n_r;
%      JacEv=JacEv + n_J;
%
%      %ResEv=ResEv + Result.ResEv;
%      %JacEv=JacEv + Result.JacEv;
%
%      Iter =Iter  + Result.Iter;
%      PrintResult(Result,PriLev);
%
%      fprintf('\n');
%      xprinte(Result.x_k,'x_k:');
%      fprintf('ResEv %d - JacEv %d - Iter %d\n',n_r,n_J,Result.Iter);
%      fprintf('Minimax L2-result %30.20f\n',max(Result.r_k));
%
%      xNew = Result.x_k;
%      diff=xNew-xOld;
%
%      if all(abs(diff) <= epsX)
%         notConv=0;
%      end
%
%      xOld=xNew;
%
%      Iter = Iter + 1;
%   end
%   fprintf('\n');
%   fprintf('Minimax L2-result %30.20f\n',max(Result.r_k));
%   fprintf('Total ResEv %d - Total JacEv %d - Total Iter %d\n',...
%            ResEv,JacEv,Iter);
%
% ---------------------------------------------------------------------------
%elseif InfType == 3 % NOT READY YET 
%   if isempty(Prob.x_0)
%      n        = max([length(Prob.x_L),length(Prob.x_U),size(Prob.A,2)]);
%      Prob.x_0 = zeros(n,1);
%   else
%      n        = length(Prob.x_0);
%   end
%
%   % Use LS formulation with nonlinear equality constraint 
%   r    = feval(Prob.FUNCS.r,Prob.x_0,Prob);
%   m    = length(r);
%   rMax = max(abs(r));
%
%   z=rMax;
%   s        = max(0,z-abs(r));
%   Prob.x_0 = [Prob.x_0;z;s];
%   Prob.N   = length(Prob.x_0);
%   Prob.x_L = [Prob.x_L;0;zeros(m,1)];
%   Prob.x_U = [Prob.x_U;rMax;rMax*ones(m,1)];
%
%   Prob.c_L = [Prob.c_L;zeros(n,1)];
%   Prob.c_U = [Prob.c_U;zeros(m,1)];
%
%   if ~isempty(Prob.A)
%      Prob.A = sparse([Prob.A,zeros(size(Prob.A,1),m+1)]);
%   end
%
%   Prob.r_k = [r-z*ones(m,1)+s;z];
%
%   Prob     = tomFiles(Prob,'ls_f','ls_g','ls_H', ...
%                     'mmxLS_c','mmxLS_dc',[], 'mmxLS_r','mmxLS_J');
%
%
%   Result=tomRun(Solver,Prob,2);
%
%   fprintf('\n');
%   xprinte(Result.x_k,'x_k:');
%
%   %fprintf('ResEv %d - JacEv %d - Iter %d\n',n_r,n_J,Result.Iter);
%   fprintf('Minimax L2-result %30.20f\n',max(Result.r_k));
%
%
% ---------------------------------------------------------------------------
%else
   % Standard approach - InfType = 1

   Prob = iniSolve(Prob,solvType,2,2);
   n = Prob.N;

   % Safe-guard starting point
   Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));
   
   % Check the residual and Jacobian supplied
   Prob.Mode   = 2;
   Prob.nState = 1;
   r  = nlp_r(Prob.x_0,Prob);
   m  = length(r);
   Prob.minimax.m   = m;
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
      Prob.A = [sparse(Prob.A),sparse(mA,1)];
   else 
      mA = 0;
   end

   % Check and adjust bounds on nonlinear constraints
   ML           = length(Prob.c_L);
   MU           = length(Prob.c_U);
   M            = max(ML,MU);
   [M1,M2]      = size(Prob.ConsPattern);
   ConsDiff     = Prob.ConsDiff;
   ADCons       = Prob.ADCons;
   if M > 0 | M1 > 0
      Prob.Mode = 2;
      c         = nlp_c(Prob.x_0,Prob);
      M         = length(c);
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
         Prob.ConsPattern = [sparse(Prob.ConsPattern),sparse(M,1); ...
              sparse(Prob.JacPattern),sparse(ones(m,1))];
      else
         m = length(Prob.LS.y);
         Prob.ConsPattern = [sparse(Prob.ConsPattern),sparse(M,1); sparse(ones(m,n+1))];
      end
   elseif ~isempty(Prob.JacPattern)
      %m = size(Prob.JacPattern,1);
      if M > 0
         Prob.ConsPattern = [sparse(ones(M,n)),sparse(M,1); ...
             sparse(Prob.JacPattern),sparse(ones(m,1))];
      else
         Prob.ConsPattern = [sparse(Prob.JacPattern),sparse(ones(m,1))];
      end
   end

   % Reformulate problem
   Prob.JacPattern = [];

   rMax = max(abs(r));
   % Add one extra variables
   Prob.x_0 = [Prob.x_0;rMax];
   Prob.x_L = [Prob.x_L;Prob.f_Low];
   if isfield(Prob,'f_Upp')
      U = Prob.f_Upp;
      if isempty(U)
         Prob.x_U = [Prob.x_U;1E300];
      else
         Prob.x_U = [Prob.x_U;max(Prob.x_L(end),U)];
      end
   else
      Prob.x_U = [Prob.x_U;Inf];
   end
   Prob.N   = n + 1;
   % Integer variables
   IntVars  = DefPar(Prob.MIP,'IntVars',[]);

   % Logical vector for integers
   IV = false(n,1);

   if isempty(IntVars)
      % No binary variables B or integer variables of type I
   elseif any(IntVars==0) | all(IntVars==1)
      % Assume binary logical vector given
      IV(1:length(IntVars)) = logical(IntVars);
   else
      if any(IntVars < 1 | IntVars > n)
         error('infSolve: Illegal IntVars vector');
      end
      IV(IntVars)=1;
   end
   Prob.MIP.IntVars = find(IV); % Generate indices, works for n+1 variables also

   % Add m extra constraints from the m residuals
   Prob.c_L     = [Prob.c_L;-Inf*ones(m,1)];
   Prob.c_U     = [Prob.c_U;zeros(m,1)];
   % Must change the number of nonlinear constraints in Prob
   Prob.mNonLin = length(Prob.c_L);

   if PriLev > 0
      fprintf('The problem has %d variables (Columns),',n);
      fprintf(' infSolve adds unbounded variable z as (Column) %d\n',Prob.N);
      fprintf('This extra variable z is the objective function value, where');
      fprintf(' %d <= z <= %d\n',Prob.x_L(end),Prob.x_U(end));
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
      fprintf('The problem has %d residuals, by infSolve ',m);
      fprintf('defined as constraint (Row) %d to %d',mA+M+1,mA+M+m);
      fprintf('\n');
   end

   Prob = tomFiles(Prob,'mmx_f','mmx_g','mmx_H','mmx_c','mmx_dc');

   Prob.nState           = 0;
   Prob.minimax.args     = args;
   Prob.minimax.NumDiff  = NumDiff;
   Prob.minimax.ConsDiff = ConsDiff;
   Prob.minimax.ADObj    = ADObj;
   Prob.minimax.ADCons   = ADCons;
   Prob.NumDiff          = 0;
   Prob.ConsDiff         = 0;
   Prob.ADObj            = 0;
   Prob.ADCons           = 0;

   Result=tomRun(Solver,Prob,PriLev);
   Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation

   % Return only original number of variables, remove f(x) value slack 
   if ~isempty(Result.x_0)
      Result.x_0     = Result.x_0(1:n);
   end
   if ~isempty(Result.xState)
      Result.xState  = Result.xState(1:n);
   end
   if ~isempty(Result.cState)
      Result.cState  = Result.cState(1:M);
   end
   if ~isempty(Result.v_k)
      if length(Result.v_k) > n
         Result.v_k  = Result.v_k([1:n,n+2:n+1+M]);
      else
         Result.v_k = [];
      end
   end
   if ~isempty(Result.x_k)
      Result.x_k  = Result.x_k(1:n,1);
   end
   if ~isempty(Result.c_k)
      Result.r_k  = Result.c_k(M+1:M+m);
      Result.c_k  = Result.c_k(1:M);
   end
   if ~isempty(Result.cJac)
      Result.J_k  = Result.cJac(M+1:M+m,1:n);
      Result.cJac = Result.cJac(1:M,1:n);
      Result.g_k  = Result.J_k'*Result.r_k;
   else
      Result.g_k  = [];
   end
   Result.H_k  = [];
%end

% MODIFICATION LOG:
%
% 990522  hkh  Written.
% 020328  hkh  Revised and modified for v3.1
% 020403  hkh  Add size checking and printing. Additional input PriLev.
% 020407  hkh  Add iniSolve call to clear global variables
% 020409  hkh  Set Prob.Mode and Prob.nState before calling nlp_r
% 020409  hkh  Return as Result.x_k only original number of variables.
% 020413  hkh  Set Prob.LargeScale = 1 if a ConsPattern or JacPattern is set
% 020413  hkh  Avoid calling feval to get rMax first, must check problem
% 020416  hkh  Fix handling of numerical derivatives, send Cons-/NumDiff
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define field mLin
% 040728  med  tomFiles used instead
% 040809  med  Pragmas added for MATLAB Compiler
% 041123  hkh  Change call to tomRun
% 041130  hkh  Prob.ConsDiff must be 0, not max(NumDiff,ConsDiff);
% 041201  hkh  Force a check in ProbCheck, setting Prob.CHECK=0
% 041222  med  Safeguard added for x_0
% 050421  hkh  Comments on all but InfType = 1, other InfType not working well
% 050421  hkh  Reset fields mLin and mNonLin, otherwise mex solvers fails
% 050421  hkh  Comments on some unnecessary tests
% 050422  hkh  Avoid Lagrange multipliers in PrintResult
% 050901  med  Error if incorrect length for c_L or c_U
% 050901  med  Removed size checks on x_L, x_U, x_0, A, b_L, b_U
% 050901  med  Replaced all applicable zeros with sparse
% 050901  med  Replaced all applicable ones with sparse(ones)
% 050901  med  Avoid recalculating size(Prob.A,1)
% 050901  med  Removed setting mLin, number of constr not changed
% 060406  med  Now handling IntVars, possible to use glcDirect
% 060406  med  v_k length checked before setting
% 060814  med  FUNCS used for callbacks instead
% 060818  hkh  Use f_Low directly in Prob.x_L = [Prob.x_L;Prob.f_Low];
% 060818  hkh  Add comments about f_Low, f_Upp (set default 1E300)
% 061211  med  Updated printing (%f switched to %g)
% 070222  hkh  Revise IntVars handling, use new format
% 080604  hkh  Safeguard NumDiff/ConsDiff, if [] or bad value
% 110722  hkh  Numeric & AD differentation for both residual and constraints
% 110724  hkh  Test ConsDiff,NumDiff for ~=6, not < 6
