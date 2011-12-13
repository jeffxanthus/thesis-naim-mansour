% TOMLAB GENO MINLP Multiobjective Genetic Solver
%
% GENO solves general constrained mixed integer uni- and multi-objective
% optimization problems using a genetic algorithm (General Evolutionary
% Numerical Optimizer)
%
% GENO solves uni-objective static optimization problems of the form:
%
%  min   f(x)
%   x
%
%  s/t   x_L <=   x  <= x_U
%        b_L <= A x  <= b_U
%        c_L <= c(x) <= c_U
%        x(i) integer, for i in I
%
% and multi-objective problems of the form:
%
%  min   J(x) = { f1(x), f2(x), f_m(x) }
%   x
%
%  s/t   x_L <=   x  <= x_U
%        b_L <= A x  <= b_U
%        c_L <= c(x) <= c_U
%        x(i) integer, for i in I
%
% Uni-objective problems is best formulated in the 'glc' format.
%
% Multi-objective problems should be given in the TOMLAB format for
% constrained least squares problems, 'cls', with each element in the residual
% vector corresponding to an objective function. See clsAssign.
%
% In either case, it is recommended to keep the variable box as small as
% possible (x_L, x_U) and especially to avoid unbounded variables if possible.
%
% Calling syntax:
%
% function Result = genoTL(Prob)
%
% INPUT PARAMETERS:
%
% Prob    Problem structure in TOMLAB format. The following fields are used:
%
%   x_0       Initial vector.
%
%   x_L, x_U  Bounds on variables.
%
%   A         Linear constraint matrix.
%   b_L, b_U  Bounds on linear constraints.
%
%   c_L, c_U  Bounds on nonlinear constraints.
%
%   PriLevOpt Print level in GENO solver. A nonzero value will produce output to
%             the MATLAB command window. Note that GENO operates on a
%             transformed set of constraints <= 0.0.
%
% ---------------------------------------
% optParam    Structure in Prob, Prob.optParam.
% ---------------------------------------
%             Defines optimization parameters. Fields used:
%
%   MaxIter   Controls the maximum number of evolutionary generations. NOTE:
%             every generation involves a large number of function evaluations.
%             It is strongly recommended to try a relatively small value here,
%             approximately 100-500.
% ---------------------------------------
% MIP         Structure in Prob, Prob.MIP
% ---------------------------------------
%             Defines integer optimization parameters. Fields used:
%   IntVars:
%             If empty, all variables are assumed non-integer
%             If islogical(IntVars) (=all elements are 0/1), then
%             1 = integer variable, 0 = continuous variable.
%             If any element >1, IntVars is the indices for integer variables
%
% ---------------------------------------
% GENO        Structure with GENO solver specific fields:
% ---------------------------------------
%
%  popsize    Defines the size of the evolutionary population. Numerical experiments
%             on a variety of problems seem to suggest that the optimal population
%             size lies in the 20 – 30 range for most problems. The default
%             value if nothing is given is 20.
%
%  PrintFile  Name of a file for printing GENO progress information. The default
%             is '', i.e. no file is created.
%
%  options    Structure in Prob.GENO, Prob.GENO.options with solver specific flags and
%             scalar parameters. The following fields may be set (default values are listed):
%
%  Real-valued options
%
%    m_rate         0.05
%    bm_rate        0.005
%    p_u_xover      0.00   Probability threshold for the uniform cross-over operator
%    p_s_xover      0.55   Probability threshold for the simple cross-over operator
%    p_a_xover      0.55   Probability threshold for the arithmetic cross-over operator
%    p_b_xover      0.00   Probability threshold for the boundary cross-over operator
%    p_h_xover      0.55   Probability threshold for the heuristic cross-over operator
%    p_d_xover      0.55   Probability threshold for the differential cross-over operator
%    d_factor       0.80   Weighting factor on the direction component of the differential cross-over operator
%   			
%    quantum_0      0.1    Specifies the initial size of quanta. In setting this parameter, the object
%                          should be to ensure that the initial population is sufficiently diverse on all
%                          dimensions. In this regard, a choice of the smaller between 0.1 and 10% of the
%                          smallest variable range is normally efficient. If an integer solution is desired,
%                          this parameter should be set to 1.
%
%    s_p            0.965  Rank fitness selection pressure
%    ss_prob        0.50   Stochastic climb operator
%    time_constant 30.0    Stochastic climb operator
%
%  Boolean (0/1) flags:
%
%    closed_loop        0  Set to one for state feedback solution
%    cross_bs           1  Set to one to enable crossing of best chromosome
%    error_tracking     0  Set to one for dynamic tracking problems
%    loop_timer         1  Set to one to enable timing in GENO
%    maximise           0  Set to one to solve maximisation problem
%    mutate_bs          1  Set to one to mutate best chromosome
%    mutate_cb          1  Set to one to mutate "current" best chromosome
%    pareto_group       0  Set to one for grouped Nash equilibrium
%    proximity          0  The proximity parameter effectively switches GENO to a "next-to-ideal" mode of
%                          operation in the quest for a multi-objective solution to the problem.
%    solution_check     1  Specifies whether the contents of the solution matrix is displayed at the end
%                          of the program run.
%    constraints_check  1  Whether or not to display values of the constraints at the end of the program run.
%    view_vars          0  Specifies whether the output is variables values (1) or objective function values (0)
%
%  Integer-valued options (Warning: do not use int32 or similar datatypes)
%
%    rand_seed   240657    Seed value for random number generator
%    		
%    partition_point  0    Partition point for Pareto groups
%    pertub_freq_m    5    Frequency of mutating ep point
%    pertub_freq_ax   5    Frequency of ax_crossing ep point
%    pertub_freq_dx   5    Frequency of dx_crossing ep point
%    pertub_freq_hx   5    Frequency of hx_crossing ep point
%    print_number    10    Print frequency: setting the value k will print progress information every k:th generation
%
%  Character-valued options:
%
%    solution_type  'e'
%
%                  'e': Nash equilibrium solution
%                  'p': Pareto-efficient solution
%                  'i': Euclidian compromise solution
%                  'm': Nash-Pareto mixed solution
%
%    adj_mode       's' Specifies GENO search procedure adjustment mode:
%
%                  's': Singly-rational adjustment mode
%                  'g'. Group-rational adjustment mode
%
%
%
%
% OUTPUT PARAMETERS:
%
% Result      Structure with results from optimization
%
%   x_k       Matrix with optimal points
%
%   f_k       For uni-objective problems: optimal function value
%             For multi-objective problems: f_k = 0.5*r_k'*r_k, see r_k below
%
%   r_k       For uni-objective problems: not set by genoTL
%             For multi-objective problems: the vector of objective functions
%
%   Ax        Linear constraints vector at x_k
%   c_k       Nonlinear constraints vector at x_k
%
%   FuncEv    Number of function evaluations
%   ConstrEv  Number of constraint evaluations
%
%   ExitFlag  0 = Normal termination, max number of iterations reached.
%

% Anders Goran, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2006-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Aug 14, 2006.  Last modified Jul 17, 2009.

function Result = genoTL(Prob)

if nargin < 1
   error('genoTL requires one input, the Prob structure')
end

Prob = ProbCheck(Prob,'glcDirect',10,10);

PriLev = Prob.PriLevOpt;

% Unless an LS type problem, set 'glc'.. need to detect this somehow, for
% multi-objective

ptyp = Prob.probType;
lstyp = [checkType('ls'),checkType('lls'),checkType('cls')];
if( any(ptyp==lstyp) )
   % This is a multi-objective problem
   typ = ptyp;
else
   typ = checkType('glc');
end

Prob.solvType = typ;

Prob = iniSolve(Prob,typ,0,0);

Prob = ConvertConstraints(Prob);

Result = ResultDef(Prob);
Result.Solver = 'GENO';
Result.SolverAlgorithm = 'Genetic Static Optimization';

[bl,bu,n,m1,m2] = defblbu(Prob);

xl = bl(1:n);
xu = bu(1:n);

order  = n;
plan   = 1;
eqns   = 0;
grad   = 0;
cons   = Prob.m;

LS = DefPar(Prob,'LS');

if isfield(LS,'t') & isfield(LS,'y')
   agents = max( [1,length(Prob.LS.t),length(Prob.LS.y)] );
elseif isfield(Prob.LS,'y')
   agents = max(1,length(Prob.LS.y));
elseif isfield(Prob.LS,'t')
   agents = max(1,length(Prob.LS.t));
else
   agents = 1;
end


GENO    = DefPar(Prob,'GENO',[]);
popsize = DefPar(GENO,'popsize',20);

maxgens = DefPar(Prob.optParam,'MaxIter',50);

dims = [ order, agents, plan, cons, eqns, grad, maxgens, popsize ];

% State bounds
usb = xu;
lsb = xl;

% Initial and final states
x_0 = DefPar(Prob,'x_0',zeros(order,1));
is = x_0;
fs = zeros(order,1);

% Control vector bounds should be so that we can reach from x_0 (is)
% to either x_L or x_U:
ucb = xu-is;
lcb = xl-is;

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
      error('genoTL: Illegal IntVars vector');
   end
   IV(IntVars)=1;
end

% Connexion matrix
ic = []; % speye(order);

PrintFile = DefPar(GENO,'PrintFile','');

global n_f n_c

n_f = 0;
n_c = 0;

options = DefPar(GENO,'options',[]);
vars    = DefPar(GENO,'vars',[]);

MaxCPU  = DefPar(Prob,'MaxCPU',1e20);

max_time = DefPar(options,'max_time',[]);
if isempty(max_time)
  options.max_time = MaxCPU;
end

% GENO mex cannot use logical input, change to double
sol_mtx = geno(dims,lcb,ucb,lsb,usb,is,fs,double(IV),ic,vars,options,...
               PriLev,PrintFile,Prob);

Result.GENO.sol_mtx = sol_mtx;

% 2nd column contains the interesting information:
sol_mtx = sol_mtx(:,2);

if(agents==1)
   % Uni-objective
   x_k = sol_mtx(order+1:2*order);
   f_k = sol_mtx(end);
   r_k = [];
else
   % Multi-objective.
   %
   % vars(agents,order) tells which variables are found in which part of sol_mtx
   % If no vars setting exists or if it's all ones, then all variables are in
   % the first agent (and identical in the remaining agents)

   if isempty(vars) | all(vars==1)
      x_k = sol_mtx(order+1:2*order);
      r_k = zeros(agents,1);
      for k=1:agents
         r_k(k) = sol_mtx(k*(2*order+1));
      end
   else
      % Use vars to find objectives and variables
      x_k = zeros(order,1);
      r_k = zeros(agents,1);
      j = 1;
      for k=1:agents
         ix = find(vars(k,:));
         j = j + length(ix);
         x_k(ix) = sol_mtx(j:j+length(ix)-1);
         r_k(k) = sol_mtx(j+length(ix));
         j = j+length(ix)-1;
      end
   end
   % As scalar objective, take f_k = 0.5*r_k'*r_k;
   f_k = 0.5*(r_k'*r_k);
end

Result.x_k = x_k;
Result.r_k = r_k;
Result.f_k = f_k;

% Restore constraint values at final point
[f_k,g_k] = nlfunc(Result.x_k,Prob);

c_k   = zeros(m2,1);
cixEQ = Prob.cixEQ;
if ~isempty(cixEQ)
   c_k(cixEQ)    = g_k(cixEQ) + Prob.c_L(cixEQ);
end
cixLow = Prob.cixLow;
if ~isempty(cixLow)
   c_k(cixLow)    = g_k(cixLow) + Prob.c_L(cixLow);
end
cixUpp = Prob.cixUpp;
if ~isempty(cixUpp)
   c_k(cixUpp)    = Prob.c_U(cixUpp) - g_k(cixUpp);
end

Result.c_k  = c_k;

if m1>0, Result.Ax  = Prob.A*x_k; else Result.Ax = [];  end

% GENO gives no information - it just runs until it stops.

Result.Inform   = 0;
Result.ExitFlag = 0;
Result.ExitText = 'Maximum number of generations reached';

Result = endSolve(Prob,Result);


% ------------------------------------------------------------------------------

% function Prob = ConvertConstraints(Prob)
%
% The Tomlab constraints on standard form:
% b_L <= A * x <= b_U
% c_L <= c(x)  <= c_U
%
% are transformed to the GENO format, which are C(x) <= 0.0
%
% The total number of constraints after conversion is Prob.m
%
% This routine is based on ksDef.m but does not treat equalities separately.
%
% The fields set to be used in Prob are:
% AixEQ, AixLow, AixUpp
% cixEQ, cixLow, cixUpp
% m, mEQ, mIN

function Prob = ConvertConstraints(Prob)

% Never identify any constraints as equalities
Prob.AixEQ = [];
Prob.cixEQ = [];

if ~isempty(Prob.A)
   if isempty(Prob.b_L)
      Prob.AixUpp = find(~isinf(Prob.b_U));
      Prob.AixLow = [];
   elseif isempty(Prob.b_U)
      Prob.AixLow = find(~isinf(Prob.b_L));
      Prob.AixUpp = [];
   else
      Prob.AixLow = find(~isinf(Prob.b_L));
      Prob.AixUpp = find(~isinf(Prob.b_U));
   end
else
   Prob.AixLow = [];
   Prob.AixUpp = [];
end

if ~(isempty(Prob.c_L) & isempty(Prob.c_U))
   if isempty(Prob.c_L)
      Prob.cixUpp = find(~isinf(Prob.c_U));
      Prob.cixLow = [];
   elseif isempty(Prob.c_U)
      Prob.cixLow = find(~isinf(Prob.c_L));
      Prob.cixUpp = [];
   else
      Prob.cixLow = find(~isinf(Prob.c_L));
      Prob.cixUpp = find(~isinf(Prob.c_U));
   end
else
   Prob.cixLow = [];
   Prob.cixUpp = [];
end

Prob.mEQ = 0; % length(Prob.AixEQ) + length(Prob.cixEQ);
Prob.mIN = length(Prob.AixLow) + length(Prob.AixUpp) + ...
           length(Prob.cixLow) + length(Prob.cixUpp);
Prob.m   = Prob.mEQ + Prob.mIN;
Prob.sense = -1;

% MODIFICATION LOG
%
% 060814 ango Wrote file
% 060825 ango First release candidate
% 060829 ango Add MaxCPU handling
% 070222 hkh  Revise IntVars handling, use new format
% 070223 hkh  Boolean input to mex must be given as double vector
% 080229 ango Fix control bounds to correctly use initial state vales
% 090717 med  f_0 calculation updated
