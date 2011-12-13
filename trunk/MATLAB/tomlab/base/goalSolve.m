% goalSolve.m:
%
% Finds a solution to a multi-objective goal attainment optimization problem,
% with linear and nonlinear constraints.
% with the use of any suitable TOMLAB solver
%
% function Result = goalSolve(Prob, PriLev)
%
% Minimization problem:
%
%        min  max {lam: r(x) - w.*lam <= g}, where r(x) is in R^m
%         x
%        s/t   x_L <=   x  <= x_U, x is in R^n
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% g are m goals, w are m weights
%
% The goal attainment problem is solved in goalSolve by rewriting the problem
% as a general constrained optimization problem.
% One additional variable z, stored as x(n+1), is added
%
%        min    z,    where r(x) is in R^m
%         x
%        s/t   x_L <=   x(1:n)         <= x_U
%             -Inf <=   z              <= Inf
%              b_L <= A x              <= b_U
%              c_L <= c(x)             <= c_U
%                g <= r(x) - z*w       <= g,  w,g is in R^m, index 1:GE
%             -Inf <= r(x) - z*w       <= g,  w,g is in R^m, index GE+1:m
%
% GE = Prob.GoalsExact, the number of goals to be exactly fulfilled
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
%              If PriLev == 0, goalSolve is silent, except for error messages.
%              If > 0, goalSolve prints summary output about problem
%              transformation
%              goalSolve calls PrintResult(Result,PriLev), i.e. printing in
%              PrintResult is made if PriLev > 0.
%              PriLev == 2 displays standard output in PrintResult.
%
%   Use Prob = probInit('name of file',problem_number'); if solving
%   a predefined problem in the Init File (IF) format.
%
%   Extra fields used in Prob:
%
%   SolverInf   Name of the TOMLAB solver. Valid names are:
%               conSolve, nlpSolve, sTrustr, clsSolve
%               If TOMLAB /SOL is installed: minos, snopt, npsol
%
%   f_Low       A lower bound on the optimal function value.
%               Not crucial, if not set default -1E300 is used
%   f_Upp       An upper bound on the optimal function value.
%               Not crucial, if not set default 1E300 is used
%
%   LS.w        Weights w, should have length m = length(r(x))
%   LS.g        Goals g, should have length m = length(r(x))
%
%   GoalsExact  The 1st "GoalsExact" goals should be fulfilled with equality
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
%   goalSolve will create the new Prob.ConsPattern to be used by the solver
%   using the information in ConsPattern and JacPattern.
%
%
% OUTPUT PARAMETERS
% Result Structure with results from optimization. See help for the used solver
%        The output in Result, i.e. fields Result.x_k, Result.r_k, Result.J_k,
%        Result.c_k, Result.cJac, Result.x_0, Result.xState, Result.cState,
%        Result.v_k, is transformed back to the original problem.
%        Result.g_k is Result.J_k'*Result.r_k.
%        The output in Result.Prob is the result after goalSolve
%        transformed the problem, i.e. the altered Prob structure

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written Sept 28, 2002.   Last modified July 24, 2011.

function Result = goalSolve(Prob, PriLev, varargin)

%#function mmx_f mmx_g mmx_H goal_c goal_dc

if nargin < 2
    PriLev = [];
    if nargin < 1
        error('goalSolve needs input structure Prob');
    end
end
if isempty(PriLev), PriLev = 2; end

if isfield(Prob,'SolverInf')
    Solver=Prob.SolverInf;
    Solver=deblank(Solver);
else
    Solver = [];
end
if ~isfield(Prob.LS,'w')
    Prob.LS.w = [];
end
if ~isfield(Prob.LS,'g')
    Prob.LS.g = [];
end

GoalsExact = DefPar(Prob,'GoalsExact',0);
Prob.GoalsExact = GoalsExact;

solvType=checkType('con');
if ~isempty(Prob.JacPattern),  Prob.LargeScale = 1; end
if ~isempty(Prob.ConsPattern), Prob.LargeScale = 1; end
if isempty(Solver), Solver=GetSolver('con',Prob.LargeScale,0); end

Prob.CHECK = 0; % Force a test in ProbCheck
Prob=ProbCheck(Prob,Solver,solvType);
Prob=iniSolve(Prob,checkType('cls'),1,1);

Prob.minimax.r  = Prob.FUNCS.r;
Prob.minimax.J  = Prob.FUNCS.J;
Prob.minimax.c  = Prob.FUNCS.c;
Prob.minimax.dc = Prob.FUNCS.dc;

args = zeros(9,1);
if ~isempty(Prob.FUNCS.r)
    args(7)     = xnargin(Prob.FUNCS.r);
end
if ~isempty(Prob.FUNCS.J)
    args(8)     = xnargin(Prob.FUNCS.J);
end
if ~isempty(Prob.FUNCS.c)
    args(4)     = xnargin(Prob.FUNCS.c);
end
if ~isempty(Prob.FUNCS.dc)
    args(5)     = xnargin(Prob.FUNCS.dc);
end

% Standard approach - InfType = 1

n = Prob.N;

% Should not call iniSolve here, no need
% Prob = iniSolve(Prob,solvType,2,2);

% Safe-guard starting point
Prob.x_0    = max(Prob.x_L,min(Prob.x_0,Prob.x_U));

% Check the residual and Jacobian supplied
Prob.Mode   = 2;
Prob.nState = 1;

r  = nlp_r(Prob.x_0,Prob,varargin{:});
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
mmw = length(Prob.LS.w);
if mmw == 0
    Prob.LS.w = ones(m,1);
elseif mmw ~= m
    fprintf('Number of residual elements is %d\n',m);
    fprintf('The number of weights in Prob.LS.w is %d\n',mmw);
    error('Should be equal, illegal size of vector Prob.LS.w');
else
    Prob.LS.w = Prob.LS.w(:);
end
mmg = length(Prob.LS.g);
if mmg == 0
    Prob.LS.g = zeros(m,1);
elseif mmg ~= m
    fprintf('Number of residual elements is %d\n',m);
    fprintf('The number of goals in Prob.LS.g is %d\n',mmg);
    error('Should be equal, illegal size of vector Prob.LS.g');
else
    Prob.LS.g = Prob.LS.g(:);
end
Prob.Mode   = 1;
NumDiff     = Prob.NumDiff;
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
ML = length(Prob.c_L);
MU = length(Prob.c_U);
M  = max(ML,MU);
[M1,M2] = size(Prob.ConsPattern);
ConsDiff     = Prob.ConsDiff;
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
        Prob.ConsPattern = sparse([Prob.ConsPattern,zeros(M,1); ...
            Prob.JacPattern,ones(m,1)]);
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
Prob.HessPattern = sparse(Prob.N,Prob.N);

% Add m extra constraints from the m residuals
if GoalsExact > 0
    Prob.c_L = [Prob.c_L;Prob.LS.g(1:GoalsExact);-Inf*ones(m-GoalsExact,1)];
else
    Prob.c_L = [Prob.c_L;-Inf*ones(m,1)];
end
Prob.c_U = [Prob.c_U;Prob.LS.g];
Prob.mNonLin = length(Prob.c_U);

if PriLev > 0
    fprintf('The problem has %d variables (Columns),',n);
    fprintf(' goalSolve adds unbounded variable z as (Column) %d\n',Prob.N);
    fprintf('This extra variable z is the objective function value, where');
    fprintf(' %e <= z <= %e\n',Prob.x_L(end),Prob.x_U(end));
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
    fprintf('The problem has %d residuals, by goalSolve ',m);
    fprintf('defined as constraint (Row) %d to %d',mA+M+1,mA+M+m);
    fprintf('\n');
end

Prob = tomFiles(Prob,'mmx_f','mmx_g','mmx_H','goal_c','goal_dc');

Prob.nState           = 0;
Prob.minimax.args     = args;
Prob.minimax.NumDiff  = NumDiff;
Prob.minimax.ConsDiff = ConsDiff;
Prob.NumDiff          = 0;
Prob.ConsDiff         = 0;

Result=tomRun(Solver,Prob,PriLev);

Result.Prob.PrintLM =0;  % Avoid Lagrange multiplier computation
% Return only original number of variables, remove f(x) value slack

if ~isempty(Result.x_0)
    Result.x_0      = Result.x_0(1:n);
end
if ~isempty(Result.xState)
    Result.xState   = Result.xState(1:n);
end
if ~isempty(Result.cState)
    Result.cState   = Result.cState(1:M);
end

Result.Prob.x_L = Prob.x_L(1:n);
Result.Prob.x_U = Prob.x_U(1:n);

if ~isempty(Result.v_k)
    Result.v_k  = Result.v_k([1:n,n+2:n+1+M]);
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

% MODIFICATION LOG:
%
% 020928  hkh  Created, based on infSolve.m
% 030127  hkh  Return correct length on Result.Prob.x_L and x_U
% 030611  ango Correct comments
% 040111  hkh  Change call to inisolve
% 040125  hkh  Define fields mLin and mNonLin
% 040215  hkh  Bug defining mNonLin
% 040421  med  empty check for glcFast as Inf solver added.
% 040526  hkh  Avoid call to iniSolve. Set HessPattern as sparse 0 matrix
% 040609  med  isempty check on Result.x_0 added.
% 040728  med  tomFiles used instead
% 040809  med  Pragmas added for MATLAB Compiler
% 041123  hkh  Change call to tomRun
% 041124  hkh  Prob.ConsDiff must be 0, not max(NumDiff,ConsDiff);
% 041124  hkh  Use the 1st n values of the 1st x_k solution, Result.x_k(1:n,1);
% 041130  hkh  Call iniSolve to clear NARG and other globals
% 041201  hkh  Force a check in ProbCheck, setting Prob.CHECK=0
% 041222  med  Safeguard added for x_0
% 050421  hkh  Comments on some unnecessary tests
% 050422  hkh  Avoid Lagrange multipliers in PrintResult
% 050901  med  Removed size checks on x_L, x_U, x_0
% 050901  med  nU < M changed to MU < M
% 050901  med  All inputs to ConsPattern now sparse
% 050901  med  Error when c_L, c_U and c not same length
% 050901  med  Removed setting mLin, number of constr not changed
% 060704  med  Spelling updated
% 060814  med  FUNCS used for callbacks instead
% 060818  hkh  Use f_Low directly in Prob.x_L = [Prob.x_L;Prob.f_Low];
% 060818  hkh  Add comments about f_Low, f_Upp (set default 1E300)
% 110724  hkh  Test ConsDiff,NumDiff for ~=6, not < 6
