% TOMLAB /CPLEX LP, QP, MILP, MIQP and MIQQ Solver
%
% cplexTL converts the problem from the TOMLAB structure format and
% calls cplex.m. On return converts the result to the TOMLAB structure
% format.
% Also see the help for cplex.m
%
% cplexTL.m solves the following mixed integer (linear or quadratic)
% programming problem (LP, QP, MILP, MIQP, MIQQ):
%
%   minimize   0.5 * x'*F*x + c'x
%      x
%
%   subject to  x_L <=    x   <= x_U
%               b_L <=   Ax   <= b_U
%
%   and also, optionally, subject to the quadratic constraints
%
%       x'*Q_i*x + a_i'*x <= r_U(i), i = 1,2,...,n_qc
%
%   and/or subject to the logical constraints
%
%       y(j) == 1 -> d_ij'*x_j <= rhs(i),  j=1,n, i=1,n_log
%
%   where
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the CPLEX sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense.
%   If F is empty, an LP or MILP problem is solved
%   Q_i are sparse n x n Matlab matrices, a_i has dimension n and r_U(i) is a
%   scalar upper bound for the i:th quadratic constraint
%
%   Some or all x may be integer valued as specified by other input
%   variables.
%
% ---------------------------------------------------------------------------
%
% function Result = cplexTL(Prob)
%
% INPUT:
%
% Prob        Problem structure in TOMLAB format.
%             Use lpAssign, mipAssign, qpAssign, miqpAssign to
%             define the Prob structure.
%
%             Fields used in input structure Prob:
%
%   x_L, x_U  Lower and upper bounds on variables, size n x 1
%   b_L, b_U  Lower and upper bounds on linear constraints, size m x 1
%   A         Linear constraint matrix, dense or sparse m x n matrix
%
%             NOTE - all bounds vectors - if [], +/- Inf is assumed
%
%   QP.c      Linear objective function coefficients, size n x 1
%   QP.F      Quadratic matrix of size n x n
%
%   QP.qc     Structure array defining quadratic constraints ("qc").
%
%             Please note that CPLEX only handles single-sided bounds
%             on qc's. An arbitrary number of qc's is set using the Prob.QP.qc
%             structure array:
%
%          qc(1).Q   = sparse( <quadratic coefficient nxn matrix> );
%          qc(1).a   = full  ( <linear coefficient nx1 vector   > );
%          qc(1).r_U = <scalar upper bound> ;
%
%             And similarly for qc(2), ... , qc(n_qc).
%
%             The standard interpretation is x'*Q*x + c'*x <= r_U, but it is
%             possible to define an alternative sense x'*Q*x + c'*x >= r_L
%             by setting qc(i).sense to a nonzero value and specifying a
%             lower bound in qc(i).r_L.
%
%             Observe that the Q matrix must be sparse, non-empty and positive
%             semi-definite for all qc's. The linear coefficient vector qc(i).a
%             may be omitted or set empty, in which case all zeros are assumed.
%
%             Likewise, if a bound r_U or r_L is empty or not present, it
%             is assumed to be 0.0. Note that this is contrary to the usual
%             Tomlab standard, where an empty or omitted bound is assumed
%             to be +/- Inf. The reason is that a single-sided constraint with
%             an infinite bound would have no meaning.
%
%   PriLevOpt Print level in cplexTL, the cplex m-file and cplexmex C-interface.
%             = 0  Silent
%             = 1  Warnings and Errors
%             = 2  Summary information
%             = 3  More detailed information
%             > 10 Pause statements, and maximal printing (debug mode)
%
% Fields used in Prob.CPLEX (Structure with CPLEX specific parameters)
%
%   LogFile   Name of file to receive the CPLEX iteration and results log.
%             If empty or not present, no log is written.
%
%   SaveFile  Name of file for saving the CPLEX problem just prior to calling the
%             CPLEX solver. If empty, nothing is written. Also see the SaveMode
%             parameter below.
%
%   SaveMode  Integer flag indicating which format to use for the save file.
%             The following values are possible:
%
%          1: SAV  Binary SAV file                 (default)
%          2: MPS  MPS format, original format
%          3: LP   CPLEX LP format, original format
%          4: RMP  MPS file with generic names
%          5: REW  MPS file with generic names
%          6: RLP  LP  file with generic names
%
%             The SAV format is a binary format suitable for submission to
%             ILOG help desk.
%
%   confgrps  Conflict groups descriptor. Set this if conflict refinement is
%             desired in the case that infeasibility is detected by CPLEX.
%
%             A conflict group consists of lists of indices describing
%             which of the following entities are part of a group:
%
%             confgrps(i).lowercol   Column (variable) lower bounds
%             confgrps(i).uppercol   Column (variable) upper bounds
%             confgrps(i).linear     Linear rows
%             confgrps(i).quad       Quadratic constraints
%             confgrps(i).sos        Special ordered sets
%             confgrps(i).logical    Logical constraints
%
%             Additionally, the group's priority value may be assigned in
%
%             confgrps(i).priority
%
%             Please refer to the TOMLAB /CPLEX User's Guide for an example of
%             Conflict Refinement.
%
%   conflictFile
%             Name of file to write conflict information to. No file is written if
%             this input parameter is empty or if no such information is available.
%
%   sa        Structure telling whether and how you want CPLEX to perform a
%             sensitivity analysis (SA). You can complete an SA on the
%             objective function, right hand side vector, lower and
%             upper bounds. The sa structure contains four
%             sub structures:
%
%                 .obj, .rhs, .xl, .xu
%
%             Each one of these contain the field:
%
%                 .index
%
%             .index contain indices to variables or constraints
%             of which to return possible value ranges.
%
%             The .index array has to be sorted, ascending.
%
%             To get an SA of objective function on the four variables 120
%             to 123 (included) and variable 19, the sa structure
%             would look like this:
%
%                 sa.obj.index = [19 120 121 122 123];
%
%             The result is returned through the output parameter 'sa'.
%
%   logcon    Structure defining logical (or "indicator") constraints.
%             This is a special type of linear constraint which is included
%             in a mixed integer problem only if an associated binary
%             variable is equal to 1.
%
%             Note that when associating a variable with a logical
%             constraint, the variable in question will be forced to become
%             a binary variable; even if it was a continuous or integer
%             variable with bounds other than 0-1.
%
%             Each element logcon(i) describes one logical constraint:
%
%              y -> row'*x <= rhs    (also == and >= possible)
%
%             The following three fields (row,var,rhs) are mandatory:
%
%             logcon(i).row   A dense or sparse row vector with the same length
%                             as the number of variables in the problem.
%
%             logcon(i).var   The index of the variable y which should control
%                             whether the constraint is "active" or not.
%                             Must be less than or equal to the number of
%                             variables in the problem.
%
%             logcon(i).rhs   The scalar value of the right hand side of the i:th
%                             logical constraint.
%
%             The following fields are optional in the description of a
%             logical constraint:
%
%             logcon(i).sense Defines the sense of the i:th logical constraint:
%
%                0 or 'lt' : implies row*x <= rhs
%                1 or 'eq' : implies row*x == rhs
%                2 or 'gt' : implies row*x >= rhs
%
%             logcon(i).comp  Complement flag. The default value 0 (empty field or
%                             left out entirely) implies that the logical constraint
%                             is active when the associated variable is equal to 1.
%                             If setting the comp field to a nonzero value, the
%                             binary variable is complemented and the
%                             constraint will become active when the variable
%                             is zero.
%
%             logcon(i).name  A string containing a name for the i:th logical
%                             constraint. This is only used if a save file is
%                             written.
%
%   BranchPrio A nonnegative vector of length n. A priority order assigns a
%              branching priority to some or all of the integer variables
%              in a model. CPLEX performs branches on variables with a
%              higher assigned priority number before variables with a lower
%              priority; variables not assigned an explicit priority value
%              by the user are treated as having a priority value of 0.
%
%   BranchDir  A vector with -1, 0, 1 entries of length n. -1 forces branching
%              towards the lower end of the integer, while 1 forces branching
%              to the higher.
%
%   Tune       Flag controlling the CPLEX Parameter Tuning feature. The
%              following values are recognized:
%
%              0 - Disable Parameter Tuning (default).
%              1 - Enable Parameter Tuning and solve problem after tuning.
%              2 - Enable Parameter Tuning and return after tuning. No solution
%                  will be generated.
%
%              The non-default CPLEX parameters found during Parameter Tuning
%              is returned in Result.MIP.cpxControl together with any settings
%              given in Prob.MIP.cpxControl. The settings given as input are
%              considered as constants during the tuning process.
%
%   Presolve   Flag controlling the CPLEX Presolve feature. The following
%              values are recognized:
%
%              0 - No special treatment of presolve
%              1 - Invoke CPLEX Presolve on problem before optimization begins and
%                  create information about the changes made.
%              2 - As 1, but also return the presolved problem, including linear
%                  constraints, objective and bounds.
%
%              The results of the presolve phase (when Presolve>0) are returned
%              in Result.CPLEX.Presolve.
%
% Fields used in Prob.MIP:
%
% See the corresponding variables in cplex.m for an explanation
%
%   MIP.IntVars
%             Defines which variables are integers, of general type I or binary B
%             Variable indices should be in the range [1,...,n].
%             IntVars is a logical vector ==> x(find(IntVars > 0)) are integers
%             IntVars is a vector of indices ==> x(IntVars) are integers
%             (if [], then no integers of type I or B are defined)
%             cplex checks which variables has x_L=0 and x_U=1, i.e. binary.
%
%   MIP.PI
%             Integer variables of type Partially Integer (PI), i.e. takes an
%             integer value up to a specified limit, and any value above that
%             limit.
%             PI must be a structure array where:
%             PI.var  Vector of variable indices in the range [1,...,n]
%             PI.lim  A vector of limit values, for each of the variables
%                     specified in PI.var, i.e. for variable i,
%                     that is the PI variable with index j in PI.var:
%                     x(i) takes integer values in [x_L(i),PI.lim(j)] and
%                     continuous values in [PI.lim(j),x_U(i)].
%
%   MIP.SC    A vector with indices for the integer variables of type
%             Semi-continuous (SC), i.e. that takes either the value 0 or a
%             real value in the range [x_L(i),x_U(i)], assuming for some j,
%             i = SC(j), where i is an variable number in the range [1,n].
%
%   MIP.SI    A vector with indices for the integer variables of type
%             Semi-integer (SI), i.e. that takes either the value 0 or
%             an integer value in the range [x_L(i),x_U(i)], assuming for some j,
%             i = SI(j), where i is an variable number in the range [1,n].
%
%   MIP.sos1  A structure defining the Special Ordered Sets of Type One (sos1).
%             Assume there are k sets of type sos1, then
%             sos1(1).var is a vector of indices for variables in sos1, set 1.
%             sos1(1).row is the row number for the reference row identifying
%                         the ordering information for the sos1 set, i.e.
%                         A(sos1(1).row,sos1(1).var) identifies this information
%             sos1(2).var is a vector of indices for variables in sos1, set 2.
%             sos1(2).row is the row number for the reference row of sos1 set 2.
%             ...
%             sos1(k).var is a vector of indices for variables in sos1, setk.
%             sos1(k).row is the row number for the reference row of sos1 set k.
%
%   MIP.sos2  A structure defining the Special Ordered Sets of Type Two (sos2).
%             Specified exactly as sos1 sets, see sos1 input variable description
%
%   MIP.basis Vector containing a CPLEX Basis. If re-solving a similar problem
%             several times, this can be set to the 'basis' output argument of an
%             earlier call to cplex.m.
%
%             The length of this vector must be equal to the sum  of the number
%             of rows (n) and columns (m).
%
%             Furthermore, please note that if cpxControl.ADVIND is set to zero,
%             the advanced basis information will not be used.
%
%             The first m elements contain row basis information, with the
%             following possible values for non-ranged rows:
%
%           0 associated slack/surplus/artificial variable nonbasic at value 0.0
%           1 associated slack/surplus/artificial variable basic
%
%             and for ranged rows (both upper and lower bounded)
%
%           0 associated slack/surplus/artificial variable nonbasic at its lower bound
%           1 associated slack/surplus/artificial variable basic
%           2 associated slack/surplus/artificial variable nonbasic at its upper bound
%
%            The last n elements, i.e. basis(m+1:m+n) contain column
%            basis information:
%
%           0 variable at lower bound
%           1 variable is basic
%           2 variable at upper bound
%           3 variable free and nonbasic
%
%   MIP.xIP  Vector with MIP starting solution, if known. NaN can be used to indicate
%            missing values. Length should be equal to number of columns in problem.
%            Values of continuous variables are ignored.
%
%   MIP.cpxControl
%
%           cpxControl Structure, see the TOMLAB /CPLEX User's Guide.
%           The user only sets the fields corresponding to the parameters
%           to be changed from its default values.
%
%           Examples:
%
%             cpxControl.STARTALG = 1 Solve root node with Primal instead of Dual simplex
%             cpxControl.SUBALG   = 4 Solve subnodes with Barrier with crossover
%
%           When using CPLEX Parameter Tuning (see Prob.CPLEX.Tune), any settings
%           in cpxControl are fixed at their given values during tuning.
%
% Fields used in Prob.optParam: (Structure with optimization parameters)
%
%   MaxIter   Limit of iterations  (if not cpxControl.ITLIM is set)
%
%   MIP.callback
%
%           Logical vector defining which callbacks to use in CPLEX
%           If the i:th entry of logical vector callback is set, the corresponding
%           callback is defined. See Tomlab /CPLEX User's Guide
%           The callback calls the m-file specified in cplex.m.
%           The user may edit the existing files, or make copies and put
%           them before the predefined files in the Matlab path.
%
% ------------------------------------------------------------------------------
%
% OUTPUTS:
%
% Result   Structure with results (see ResultDef.m):
%
%   f_k      Function value at optimum
%   x_k      Solution vector. If MIP solution pool is active, multiple solutions
%            are returned in x_k with one solution per column. The first column
%            is always the best solution.
%   x_0      Initial  solution vector not known, set as empty
%   g_k      Exact gradient computed at optimum, computed as c or c + Fx
%
%   xState   State of variables.   Free==0; On lower == 1; On upper == 2; Fixed == 3;
%   bState   State of constraints. Free==0; On lower == 1; On upper == 2; Equality == 3;
%
%   v_k      Lagrangian multipliers (for bounds + dual solution vector)
%            v_k = [rc;v]. rc n-vector of reduced costs. v holds m dual variables
%
%   rc       Reduced costs. If ninf=0, last m == -v_k
%
%   ExitFlag Exit status, TOMLAB standard
%
%   Inform   CPLEX information parameter, see TOMLAB /CPLEX User's Guide
%
%      LP/QP Inform values, see help cplexStatus
%
%   Iter     Number of iterations / nodes visited
%
%   FuncEv   Number of function evaluations. Set to Iter.
%
%   GradEv   Number of gradient evaluations. Set to Iter if
%            QP/MIQP, otherwise 0.
%            FuncEv and ConstrEv set to Iter. GradEv=0.
%
%   ConstrEv Number of constraint evaluations. Set to 0.
%
%   QP.B     Basis vector in TOMLAB QP standard
%
%   Solver           Name of the solver  (CPLEX)
%   SolverAlgorithm  Description of the solver
%
% -----------------------------------------------
% Output fields in Result.MIP:
% -----------------------------------------------
%   MIP.slack     Slack variables (m x 1 vector)
%   MIP.ninf      Number of infeasibilities
%   MIP.sinf      Sum of infeasibilities
%   MIP.lpiter    Number of LP iterations
%   MIP.glnodes   Number of nodes visited
%   MIP.basis     basis status of constraints + variables, (m + n x 1 vector)
%                 in the CPLEX format, fields xState and bState has the same
%                 information in the TOMLAB format.
%   MIP.cpxControl  Structure with non-default CPLEX parameters generated during
%             CPLEX Parameter Tuning.
%
% Also set into the Result.MIP output structure is:
%
%   MIP.cpxRetVec  A vector with information on return from CPLEX, see the
%                  TOMLAB /CPLEX User's Guide for a description of each element
%
% -----------------------------------------------
% Output fields in Result.CPLEX:
% -----------------------------------------------
%
%   confstat Structure with extended conflict status information. This
%            output is a replica of the Prob.CPLEX.confgrps input argument
%            with the added fields 'status' and 'istat'. confstat(k).status
%            gives a text description of the status of conflict group k;
%            the corresponding istat field is the numeric value also
%            available in iconfstat(k).
%
%   iconfstat Conflict status information. For an infeasible problem where
%            at least one conflict group have been supplied in the confgrps
%            input argument, this output argument contains the status of
%            each conflict group, in the same order as given in the confgrps input.
%
%            The following values are possible:
%
%           -1  Excluded
%            0  Possible member
%            1  Possible member LB
%            2  Possible member UB
%            3  Member
%            4  Upper bound
%            5  Lower bound
%
%            If confstat is empty even though Conflict Refinement has been
%            requested, there was a problem in the refinement process.
%
%   sa      Structure with information about the requested SA, if requested.
%           The fields:
%
%               obj         Ranges for the variables in the objective function.
%
%               rhs         Ranges for the right hand side values.
%
%               xl          Ranges for the lower bound values.
%
%               xu          Ranges for the upper bound values.
%
%           These fields are structures themselves. All four structures
%           have identical field names:
%
%               status      Status of the SA operation. Possible values:
%
%                            1  Successful
%                            0  SA not requested.
%                           -1  Error: begin is greater than end.
%                           -2  Error: The selected range (begin...end) stretches
%                               out of available variables or constraints.
%                           -3  Error: No SA available.
%
%               lower       The lower range.
%
%               upper       The upper range.
%
%   Presolve Structure with information about the changes CPLEX Presolve made to
%            the problem before solving. The amount of information depends on
%            the value of the Prob.CPLEX.Presolve flag.
%
%            For Prob.CPLEX.Presolve=0, this output is empty.
%
%            For Prob.CPLEX.Presolve=1, the Presolve output contains six arrays,
%            pcstat, prstat, ocstat, orstat, status and objoffset.
%
%               status  Flag telling status of presolve results:
%
%                       0  Problem was not presolved or no reductions were made
%                       1  A presolved problem exists
%                       2  The original problem was reduced to an empty problem
%
%            The remaining fields of Result.CPLEX.Presolve will contain useful
%            information only if status==1.
%
%               pcstat  Contains information about variables in the original
%                       problem. For each element pcstat(i):
%
%                       >=1 variable i corresponds to variable pcstat(i) in
%                       the presolved problem
%                       -1 variable i is fixed to its lower bound
%                       -2 variable i is fixed to its upper bound
%                       -3 variable i is fixed to some other value
%                       -4 variable i is aggregated out
%                       -5 variable i is deleted or merged for some other reason
%
%               prstat  Contains information about constraints (rows) in the original
%                       problem. For each element prstat(i):
%
%                       >=1 row i corresponds to row prstat(i) in the
%                       presolved problem
%                       -1 row i is redundant
%                       -2 row i is used for aggregation
%                       -3 row i is deleted for some other reason
%
%               ocstat  Contains information about variables in the
%                       presolved problem. For each element ocstat(i):
%
%                       >=1 variable i in the presolved problem corresponds
%                       to variable ocstat(i) in the original problem.
%                       -1 variable i corresponds to a linear combination
%                       of some variables in the original problem.
%
%               orstat  Contains information about constraints (rows) in
%                       the presolved problem. For each element orstat(i):
%
%                       >=1 row i in the presolved problem corresponds to
%                       row orstat(i) in the original problem
%                       -1 row i is created by, for example, merging two
%                       rows in the original problem.

% The return vector information is also available as a global variable
% after the run, do: global cpxRetVec

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 11.2.0$
% Written July 8, 1999.      Last modified Nov 3, 2008.

function Result = cplex11TL(Prob)

%#function lp_f lp_g lp_H cpx2cbinfo cpx2retvec

Prob.MIP.cpxControl.PREIND=0;
Prob.PriLevOpt=2;

if nargin < 1
   error('cplexTL needs the Prob structure as input');
end

% Information on return from CPLEX
global cpxRetVec

global MAX_x MAX_c % Max number of variables and constraints to print

Prob.solvType = 11; % MIQP solver

Prob = iniSolveMini(Prob);

nCALLBACKS=15;
%DEBUG = 0;

PriLev=Prob.PriLevOpt;

Result=ResultDef(Prob);
Result.Solver='CPLEX';

% Initial checks on the inputs
%
% Define lower and upper bound arrays for CPLEX
%
% Inf are changed to BIG (default = 1E12), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints

% The cplex MEX has 1E12 as default, but is changed if BIG is noempty
Prob.BIG=DefPar(Prob,'BIG',1E12);

[bl, bu, n, m] = defblbu(Prob, Prob.BIG);

nTot=n+m;

% Initial values (cplex does not use x_0)

Fzero = isempty(Prob.QP.F);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   if nA~=n, error('Linear constraints A MUST have n columns!'); end
   %   fprintf('mA = %i, m = %i\n', mA, m);
   if mA~=m, error('Linear constraints A MUST have m rows!'); end
end

Result.f_0=0;

% Check if any linear part
c = Prob.QP.c(:);

if isempty(c), c=zeros(n,1); end


Prob.MIP.IntVars=DefPar(Prob.MIP,'IntVars',[]);
Prob.MIP.sos1=DefPar(Prob.MIP,'sos1',[]);
Prob.MIP.sos2=DefPar(Prob.MIP,'sos2',[]);
Prob.MIP.SC=DefPar(Prob.MIP,'SC',[]);
Prob.MIP.SI=DefPar(Prob.MIP,'SI',[]);
Prob.MIP.PI=DefPar(Prob.MIP,'PI',[]);

if isempty(Prob.MIP)
   MIP=0;

elseif isempty(Prob.MIP.IntVars)
   if isempty(Prob.MIP.sos1) & isempty(Prob.MIP.sos2) & ...
         isempty(Prob.MIP.SC) & isempty(Prob.MIP.SI)  & isempty(Prob.MIP.PI)
      MIP=0;
   else
      MIP=1;
   end
else
   MIP=1;
end

% CPLEX settings structure
cpxControl = DefPar(Prob.MIP,'cpxControl',[]);

% Generate a SolverAlgorithm text based on problem type and user setting
ptyp = checkType(Prob.probType);
switch(ptyp)
   case 'lp',
      alg = DefPar(cpxControl,'LPMETHOD',0);
      switch(alg)
         case {0,2}
            algtxt = 'CPLEX Dual Simplex LP solver';
         case 1,
            algtxt = 'CPLEX Primal Simplex LP solver';
         case 3,
            algtxt = 'CPLEX Network LP solver';
         case 4,
            algtxt = 'CPLEX Barrier LP solver';
         case 5,
            algtxt = 'CPLEX Sifting LP solver';
         case 6,
            algtxt = 'CPLEX Concurrent LP solver';
         otherwise,
            algtxt = 'CPLEX LP/QP/MILP/MIQP/MIQQ solver';
      end
      
   case 'qp',
      alg = DefPar(cpxControl,'QPMETHOD',0);
      switch(alg)
         case {0,4}
            algtxt = 'CPLEX Barrier QP solver';
         case 1,
            algtxt = 'CPLEX Primal Simplex QP solver';
         case 2,
            algtxt = 'CPLEX Dual Simplex QP solver';
         case 3,
            algtxt = 'CPLEX Network QP solver';
         case 5,
            algtxt = 'CPLEX Sifting QP solver';
         case 6,
            algtxt = 'CPLEX Concurrent QP solver';
         otherwise,
            algtxt = 'CPLEX LP/QP/MILP/MIQP/MIQQ solver';
      end
      
   case {'mip','miqp','miqq'}
      algtxt = ['CPLEX Branch-and-Cut ' upper(ptyp) ' solver' ];
      
   otherwise,
      algtxt = 'CPLEX LP/QP/MILP/MIQP/MIQQ solver';
end

Result.SolverAlgorithm = algtxt;

% Callbacks vector
callback   = DefPar(Prob.MIP,'callback',zeros(nCALLBACKS,1));

if length(callback) < nCALLBACKS
   callback=[callback(:);zeros(nCALLBACKS-length(callback),1)];
end

if MIP
   IntVars = DefPar(Prob.MIP,'IntVars',[]);
   PI      = DefPar(Prob.MIP,'PI',[]);
   SC      = DefPar(Prob.MIP,'SC',[]);
   SI      = DefPar(Prob.MIP,'SI',[]);
   sos1    = DefPar(Prob.MIP,'sos1',[]);
   sos2    = DefPar(Prob.MIP,'sos2',[]);
   xIP     = DefPar(Prob.MIP,'xIP',[]);
else
   IntVars = [];
   PI      = [];
   SC      = [];
   SI      = [];
   sos1    = [];
   sos2    = [];
   xIP     = [];
end

%Simplex iteration limit
if ~isfield(Prob.optParam,'MaxIter')
   Prob.optParam.MaxIter=2000;
end
if ~isfield(cpxControl,'ITLIM') & Prob.optParam.MaxIter~=2000
   cpxControl.ITLIM=Prob.optParam.MaxIter;
end

CPLEX = DefPar(Prob,'CPLEX');

SaveFile = DefPar(CPLEX,'SaveFile','');
SaveMode = DefPar(CPLEX,'SaveMode',1);
LogFile  = DefPar(CPLEX,'LogFile','');

if ~isempty(SaveFile) & (SaveMode<1 | SaveMode>6)
   error('SaveFile is nonempty and SaveMode is not in 1,2,...,6')
end

% Advanced starting basis, if available
advbas = DefPar(Prob.MIP,'basis',[]);

% Switch ADVIND on, if not explicitly set by the user
if length(advbas) == nTot
   ADVIND = DefPar(cpxControl,'ADVIND',[]);
   if isempty(ADVIND)
      cpxControl.ADVIND = 1;
   end
end

saRequest  = DefPar(CPLEX,'sa',[]);

% Quadratic constraints
qc = DefPar(Prob.QP,'qc',[]);

% Conflict groups information
confgrps = DefPar(CPLEX,'confgrps',[]);
conffile = DefPar(CPLEX,'conflictFile',[]);

% Logical constraint information
logcon = DefPar(CPLEX,'logcon',[]);

% Branching
branchprio = DefPar(CPLEX,'BranchPrio',[]);
branchdir  = DefPar(CPLEX,'BranchDir',[]);

% Parameter Tuning
Tune     = DefPar(CPLEX,'Tune',0);
Presolve = DefPar(CPLEX,'Presolve',0);

cpxSettings = struct('tune',Tune,'presolve',Presolve);

[x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes, ...
   confstat, iconfstat, sa, cpxControlOut, presolve ] = ...
   cplex11(c, Prob.A, bl(1:n), bu(1:n), bl(n+1:nTot), bu(n+1:nTot), cpxControl,...
   callback, PriLev, Prob, IntVars, PI, SC, SI, sos1, sos2, Prob.QP.F,...
   LogFile, SaveFile, SaveMode, qc, confgrps, conffile, saRequest, advbas, xIP, ...
   logcon, branchprio, branchdir, cpxSettings );

[Result.ExitText,Result.ExitFlag] = cplexStatus(Inform);

Result.MIP.cpxControl = cpxControlOut;

% Result printing
if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nTOMLAB / CPLEX solving Problem %d: ',Prob.P);
   fprintf('Inform %d = ',Inform);

   fprintf('%s\n',Result.ExitText);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n');

   fprintf('\nObjective function at x (obj) %25.16f\n\n',f_k);
   if MIP
      fprintf('LP iterations%7d. ',lpiter);
   else
      fprintf('Nodes visited%7d. ',glnodes);
   end
   fprintf('\n');
   fprintf('Number of infeasibilities%7d. ',ninf);
   fprintf('Sum of infeasibilities %e. ',sinf);
   fprintf('\n');
end


if PriLev > 1
   if isempty(MAX_x)
      MAX_x=length(x);
   end
   fprintf('Optimal x = \n');
   xprinte(x(1:min(length(x),MAX_x)),'x:  ');
end


if PriLev > 3
   if isempty(MAX_c)
      MAX_c=20;
   end
   fprintf('Dual variables (Lagrangian multipliers) v_k = \n');
   xprinte(v(1:min(length(v),MAX_c)),'v_k:');

   fprintf('Reduced costs rc: Last %d elements should be -v_k\n',length(v));
   xprint(rc(1:min(length(rc),MAX_c+MAX_x)),'rc: ',' %14.9f',5);
end

Result.f_k=f_k;
Result.x_0=[];
Result.x_k=x;
Result.v_k=[rc;v];

% Multiple solutions? Do the rest of the steps on the first one. 
if size(x,2)>1    
   x = x(:,1);
end

if ~isempty(c)
   if ~Fzero
      Result.g_k=Prob.QP.F*x+c;
   else
      Result.g_k=c;
   end
elseif Fzero
   Result.g_k=[];
else
   Result.g_k=Prob.QP.F*x;
end

Result.c_k=[];
Result.cJac=[];

% State variable: Free==0; On lower == 1; On upper == 2; Fixed == 3;
if ~isempty(Prob.A)
   Result = StateDef(Result, x(1:n), Prob.A(1:m,:)*x(1:n), [], ...
      Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu);
else
   Result = StateDef(Result, x(1:n), [], [], ...
      Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu);
end

if PriLev > 2
   fprintf('State vector xState for x = \n');
   xprinti(basis(1:min(length(basis),MAX_x)),'xS: ');
end
if MIP
   Result.Iter     = glnodes; % Number of nodes visited
else
   Result.Iter     = lpiter;  % Simplex iterations
end
Result.FuncEv   = lpiter;
if ~Fzero
   Result.GradEv   = lpiter;
else
   Result.GradEv   = 0;
end

Result.ConstrEv = lpiter;

if 0
   Result.Inform   = glstatus;
else
   Result.Inform   = Inform;
end

Result.MIP.ninf=ninf;
Result.MIP.sinf=sinf;
Result.MIP.slack=slack;
Result.MIP.lpiter=lpiter;
Result.MIP.glnodes=glnodes;
Result.MIP.basis=basis;

Result.CPLEX.sa = sa;
Result.CPLEX.confstat = confstat;
Result.CPLEX.Presolve = presolve;

Result.CPLEX.iconfstat = zeros(length(confstat),1);
for k=1:length(confstat)
   Result.CPLEX.iconfstat(k) = confstat(k).istat;
end

% Return information from CPLEX
Result.MIP.cpxRetVec = cpxRetVec;
if Prob.mNonLin > 0
   Result.c_k = qp_c(Result.x_k, Prob); 
end

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 020818 hkh  Last revision for xpress
% 020805 fhe  Modified xpressTL to be used with CPLEX as cplexTL
% 020922 hkh  Revision, using cpxControl for CPLEX parameters
% 030731 ango ExitFlag now conforms with TOMLAB standard, comments revision
% 030930 ango Added LogFile and SaveFile features
% 031014 ango Better check on SaveMode and SaveFile values, also commented.
% 031125 ango Adapt status messages and codes to CPLEX 9.0
% 040103 hkh  Revision for v4.2, call iniSolve and endSolve
% 040113 ango Comments revised for quadratic constraints
% 040528 hkh  Avoid lowering ITLIM, unless user has set optParam.MaxIter
% 040528 hkh  Comment out all DEBUG parts
% 040803 med  Added pragmas for MATLAB Compiler
% 040818 fhe  Added iis and sa and some help text for iis and sa.
% 040825 ango Slight change for sa help text.
% 041213 hkh  Use BIG as Prob.BIG is nonempty, otherwise 1E12
% 041213 hkh  Report glnodes as Result.Iter for MIP problems
% 041213 hkh  Use lpiter for Result.FuncEv, GradEv, ConstrEv
% 050201 ango Add support for advanced basis
% 050202 med  Removed old code and updated help
% 050414 ango mipstart added
% 050504 ango mipstart changed, use xIP instead
% 050919 med  xState, bState modified
% 060124 ango CPLEX 10 adaptations
% 060203 ango Final CPLEX 10 changes, IIS removed
% 060206 med  StateDef call modified
% 060803 ango Help updated
% 061120 ango SolverAlgorithm text improved
% 070221 hkh  IntVars format revised
% 070223 hkh  Use new routine cplexStatus to compute ExitText and ExitFlag
% 070526 med  mlint check
% 070601 ango 14 callbacks with new Incumbent Callback
% 070802 ango 15 callbacks with new User Cut Callback
% 071116 med  BranchPrio and BranchDir added
% 080201 ango Presolve feature added
% 080207 med  Help corrected