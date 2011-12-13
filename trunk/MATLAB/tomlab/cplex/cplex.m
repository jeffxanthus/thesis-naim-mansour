% TOMLAB /CPLEX LP, QP, MILP, MIQP and MIQQ Solver
%
% -----------------------------------------------------------------------
%
%   cplex.m solves the following
%   mixed integer (linear or quadratic) programming problem (MILP, MIQP):
%
%   minimize   0.5 * x'*F*x + c'x     subject to:
%      x             x_L <=    x   <= x_U
%                    b_L <=   Ax   <= b_U
%   where
%
%   A is an m x n dense or sparse Matlab matrix (linear constraints)
%   A is transformed to the CPLEX sparse matrix format.
%   c, x_L, x_U has dimension n
%   b_L, b_U has dimension m
%   F is a n x n symmetric matrix, sparse or dense.
%   If F is empty, an LP or MILP problem is solved
%
%   CPLEX also handles convex quadratic constraints (QC).
%   QC's are defined as
%
%               x'*Qj*x + aj'*x <R> rj   , j=1,2,...,nq
%
%   where <R> is one and only one of the relations <= or >=.
%
%   The aj vector may be zero for any given QC, but the Qj matrix must have at
%   least one nonzero element.
%
% ------------------------------------------------------------------------
%
% function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
%          glnodes, confstat, iconfstat, sa, cpxControl, presolve] = ...
%          cplex(c, A, x_L, x_U, b_L, b_U, ...
%          cpxControl, callback, PriLev, Prob, IntVars, PI, SC, SI, ...
%          sos1, sos2, F, logfile, savefile, savemode, qc, ...
%          confgrps, conflictFile, saRequest, basis, xIP, logcon, branchprio, ...
%          branchdir, cpxSettings);
%
% INPUT:
% c         Linear objective function cost coeffs, n x 1.
% A         Linear constraint matrix, dense or sparse m x n matrix.
%           cplex.m converts the matrix to a sparse format.
% x_L       Lower bounds on x. (if [], then x_L=0 assumed)
% x_U       Upper bounds on x
% b_L       Lower bounds on linear constraints
%
%    The following parameters are optional:
%
% b_U        Upper bounds on linear constraints (if [], then b_U=b_L assumed)
% cpxControl Structure, where the fields are set to the CPLEX
%            parameters that the
%            user wants to specify values for. The prefix CPX_PARAM
%            is not used. For example:
%
%  cpxControl.STARTALG = 1 Solve root node with Primal instead of Dual simplex
%  cpxControl.SUBALG   = 4 Solve subnodes with Barrier with crossover
%
% callback  Logical vector defining which callbacks to use in CPLEX
%  If the i:th entry of logical vector callback is set, the corresponding
%  callback is defined. See CPLEX manuals and the Tomlab /CPLEX User's Guide
%  The callback calls the m-file specified below. The user may edit this file,
%  or make a new copy, which is put before in the Matlab path.
%
% callback(1)  cpxcb_PRIM.m       Primal simplex callback
%         (2)  cpxcb_DUAL.m       Dual simplex callback
%         (3)  cpxcb_PRIMCROSS.m  Primal crossover callback
%         (4)  cpxcb_DUALCROSS.m  Dual crossover callback
%         (5)  cpxcb_BARRIER.m    Barrier log callback
%         (6)  cpxcb_PRESOLVE.m   Presolve callback
%         (7)  cpxcb_MIP.m        MIP callback
%         (8)  cpxcb_MIPPROBE.m   MIP probe or clique merging callback
%         (9)  cpxcb_FRACCUT.m    Gomory fractional cut callback
%         (10) cpxcb_DISJCUT.m    Disjunctive cut callback
%         (11) cpxcb_FLOWMIR.m    Mixed integer rounding cut callback
%         (12) cpxcb_QPBARRIER.m  QP Barrier callback
%         (13) cpxcb_QPSIMPLEX.m  QP Simplex callback
%         (14) cpxcb_INCUMBENT.m  MIP Incumbent callback
%         (15) cpxcb_USERCUT.m    MIP User cut callback
%
% PriLev    Printing level in the CPLEX m-file and CPLEX C-interface.
%           = 0    Silent
%           = 1    Warnings and Errors
%           = 2    Summary information
%           = 3    More detailed information
%
%           > 10   Pause statements, and maximal printing (debug mode)
%
% Prob      A structure. If TOMLAB calls CPLEX, then Prob is the standard
%           TOMLAB problem structure, otherwise the user optionally may set:
%           Prob.P = ProblemNumber, where ProblemNumber is some integer.
%           If any callback is defined (see description of callback) then
%           problem arrays are set as fields in Prob, and the Prob structure
%           is always passed to the callback routines as the last parameter.
%           The defined fields are Prob.QP.c, Prob.QP.F,
%           Prob.x_L, Prob.x_U, Prob.A, Prob.b_L, Prob.b_U.
%           (if input is [], then Prob.P=1 is set)
%
%           NOTE: The cplex MEX is using values >= 1E12 to define infinite
%           values in lower/upper bounds of variables and constraints
%           By setting Prob.BIG to a nonempty positive big value, this value
%           will be used by the MEX.
%           DO NOT USE the MATLAB inf or NaN value in any arrays.
%           Convert any inf value to 1E12, -inf to -1E12 (or the value set
%           as Prob.BIG).
%
% IntVars   Defines which variables are integers, of general type I or binary B
%           Variable indices should be in the range [1,...,n].
%           IntVars is a logical vector ==> x(find(IntVars > 0)) are integers
%           IntVars is a vector of indices ==> x(IntVars) are integers
%           (if [], then no integers of type I or B are defined)
%           CPLEX checks which variables has x_L=0 and x_U=1, i.e. binary.
%
% PI        Integer variables of type Partially Integer (PI), i.e. takes an
%           integer value up to a specified limit, and any value above that
%           limit.
%           PI must be a structure array where:
%           PI.var  Vector of variable indices in the range [1,...,n]
%           PI.lim  A vector of limit values, for each of the variables
%                   specified in PI.var, i.e. for variable i,
%                   that is the PI variable with index j in PI.var:
%                   x(i) takes integer values in [x_L(i),PI.lim(j)] and
%                   continuous values in [PI.lim(j),x_U(i)].
%
% SC        A vector with indices for the integer variables of type
%           Semi-continuous (SC), i.e. that takes either the value 0 or a
%           real value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SC(j), where i is an variable number in the range [1,n].
%
% SI        A vector with indices for the integer variables of type
%           Semi-integer (SI), i.e. that takes either the value 0 or
%           an integer value in the range [x_L(i),x_U(i)], assuming for some j,
%           i = SI(j), where i is an variable number in the range [1,n].
%
% sos1      A structure defining the Special Ordered Sets of Type One (sos1).
%           Assume there are k sets of type sos1, then
%           sos1(1).var is a vector of indices for variables in sos1, set 1.
%           sos1(1).row is the row number for the reference row identifying
%                       the ordering information for the sos1 set, i.e.
%                       A(sos1(1).row,sos1(1).var) identifies this information
%           sos1(2).var is a vector of indices for variables in sos1, set 2.
%           sos1(2).row is the row number for the reference row of sos1 set 2.
%           ...
%           sos1(k).var is a vector of indices for variables in sos1, setk.
%           sos1(2).row is the row number for the reference row of sos1 set k.
%
% sos2      A structure defining the Special Ordered Sets of Type Two (sos2).
%           Specified exactly as sos1 sets, see sos1 input variable description
%
% F         Square dense or sparse matrix. Empty if non-quadratic problem.
%
% logfile   Name of file to write CPLEX log to. If empty, no log is
%           written.
%
% savefile  Name of file to write CPLEX problem just prior to calling the
%           CPLEX solver. If empty, nothing is written. Also see the savemode
%           parameter below.
%
% savemode  Integer flag indicating which format to use for the save file.
%           The following values are possible:
%
%        1: SAV  Binary SAV file
%        2: MPS  MPS format, original format
%        3: LP   CPLEX LP format, original format
%        4: RMP  MPS file with generic names
%        5: REW  MPS file with generic names
%        6: RLP  LP  file with generic names
%
%            The SAV format is a binary format suitable for submission to
%            ILOG help desk.
%
% qc       Structure array defining quadratic constraints ("qc").
%
%           Please note that CPLEX only handles single-sided bounds
%           on qc's. An arbitrary number of qc's is set using the Prob.QP.qc
%           structure array:
%
%          qc(1).Q   = sparse( <quadratic coefficient nxn matrix> );
%          qc(1).a   = full  ( <linear coefficient nx1 vector   > );
%          qc(1).r_U = <scalar upper bound> ;
%
%           And similarly for qc(2), ... , qc(n_qc).
%
%           The standard interpretation is x'*Q*x + c'*x <= r_U, but it is
%           possible to define an alternative sense x'*Q*x + c'*x >= r_L
%           by setting qc(i).sense to a nonzero value and specifying a
%           lower bound in qc(i).r_L.
%
%           Observe that the Q matrix must be sparse, non-empty and positive
%           semi-definite for all qc's. The linear coefficient vector qc(i).a
%           may be omitted or set empty, in which case all zeros are assumed.
%
%           Likewise, if a bound r_U or r_L is empty or not present, it
%           is assumed to be 0.0. Note that this is contrary to the usual
%           Tomlab standard, where an empty or omitted bound is assumed
%           to be +/- Inf. The reason is that a single-sided constraint with
%           an infinite bound would have no meaning.
%
% confgrps     Conflict groups descriptor. Set this if conflict refinement is
%              desired in the case that infeasibility is detected by CPLEX.
%
%              A conflict group consists of lists of indices describing
%              which of the following entities are part of a group:
%
%              confgrps(i).lowercol   Column (variable) lower bounds
%              confgrps(i).uppercol   Column (variable) upper bounds
%              confgrps(i).linear     Linear rows
%              confgrps(i).quad       Quadratic constraints
%              confgrps(i).sos        Special ordered sets
%              confgrps(i).indicator  Logical constraints
%
%              Additionally, the group's priority value may be assigned in
%
%              confgrps(i).priority
%
%              Please refer to the TOMLAB /CPLEX User's Guide for an example of
%              Conflict Refinement.
%
% conflictFile  Name of file to write conflict information to. No file is written if
%               this input parameter is empty or if no such information is available.
%
% saRequest     Structure telling whether and how you want CPLEX to perform a
%               sensitivity analysis (SA). You can complete an SA on the
%               objective function, right hand side vector, lower and
%               upper bounds. The saRequest structure contains four
%               sub structures:
%
%                   .obj, .rhs, .xl, .xu
%
%               Each one of these contain the field:
%
%                   .index
%
%               .index contain indices to variables or constraints
%               of which to return possible value ranges.
%
%               The .index array has to be sorted, ascending.
%
%               To get an SA of objective function on the four variables 120
%               to 123 (included) and variable 19, the saRequest structure
%               would look like this:
%
%                   saRequest.obj.index = [19 120 121 122 123];
%
%               The result is returned through the output parameter 'sa'.
%
% basis         Vector with CPLEX starting basis.
%               If re-solving a similar problem several times, this can be
%               set to the 'basis' output argument of an earlier call
%               to cplex.m. The length of this vector must be equal to the sum
%               of the number of rows (m) and columns (n).
%
%               The first m elements contain row basis information, with the
%               following possible values for non-ranged rows:
%
%             0 associated slack/surplus/artificial variable nonbasic at value 0.0
%             1 associated slack/surplus/artificial variable basic
%
%               and for ranged rows (both upper and lower bounded)
%
%             0 associated slack/surplus/artificial variable nonbasic at its lower bound
%             1 associated slack/surplus/artificial variable basic
%             2 associated slack/surplus/artificial variable nonbasic at its upper bound
%
%
%               The last n elements, i.e. basis(m+1:m+n) contain column
%               basis information:
%
%             0 variable at lower bound
%             1 variable is basic
%             2 variable at upper bound
%             3 variable free and nonbasic
%
% xIP          Vector with MIP starting solution, if known. Missing values
%              may be set to NaN. Length should be equal to number of columns in problem.
%              For xIP to be taken into account by CPLEX, the parameter cpxControl.ADVIND
%              must be set to 1 or 2. Unless the caller has explicitly specified 
%              a value for ADVIND in cpxControl, cplex.m will set ADVIND=1 whenever 
%              the xIP is present and nonempty. 
%
% logcon       Structure defining logical (or "indicator") constraints.
%              This is a special type of linear constraint which is included
%              in a mixed integer problem only if an associated binary
%              variable is equal to 1.
%
%              Note that when associating a variable with a logical
%              constraint, the variable in question will be forced to become
%              a binary variable; even if it was a continuous or integer
%              variable with bounds other than 0-1.
%
%              Each element logcon(i) describes one logical constraint:
%
%               y -> row'*x <= rhs    (also == and >= possible)
%
%              The following three fields (row,var,rhs) are mandatory:
%
%              logcon(i).row   A dense or sparse row vector with the same length
%                              as the number of variables in the problem.
%
%              logcon(i).var   The index of the variable y which should control
%                              whether the constraint is "active" or not.
%                              Must be less than or equal to the number of
%                              variables in the problem.
%
%              logcon(i).rhs   The scalar value of the right hand side of the i:th
%                              logical constraint.
%
%              The following fields are optional in the description of a
%              logical constraint:
%
%              logcon(i).sense Defines the sense of the i:th logical constraint:
%
%                 0 or 'lt' : implies row*x <= rhs
%                 1 or 'eq' : implies row*x == rhs
%                 2 or 'gt' : implies row*x >= rhs
%
%              logcon(i).comp  Complement flag. The default value 0 (empty field or
%                              left out entirely) implies that the logical constraint
%                              is active when the associated variable is equal to 1.
%                              If setting the comp field to a nonzero value, the
%                              binary variable is complemented and the
%                              constraint will become active when the variable
%                              is zero.
%
%              logcon(i).name  A string containing a name for the i:th logical
%                              constraint. This is only used if a save file is
%                              written.
%
% branchprio   A nonnegative vector of length n. A priority order assigns a
%              branching priority to some or all of the integer variables
%              in a model. CPLEX performs branches on variables with a
%              higher assigned priority number before variables with a lower
%              priority; variables not assigned an explicit priority value
%              by the user are treated as having a priority value of 0.
%
% branchdir    A vector with -1, 0, 1 entries of length n. -1 forces branching
%              towards the lower end of the integer, while 1 forces branching
%              to the higher. 0 chooses the direction set by the
%              BRDIR parameter (see cpxControl)
%
% cpxSettings  Structure with flags controlling various CPLEX features.
%              Currently the following fields are defined:
%
%   tune       Flag controlling the CPLEX Parameter Tuning feature. The
%              following values are recognized:
%
%              0 - Disable Parameter Tuning (default).
%              1 - Enable Parameter Tuning and solve problem after tuning.
%              2 - Enable Parameter Tuning and return after tuning. No solution
%                  will be generated.
%
%              The non-default CPLEX parameters found during Parameter Tuning
%              is returned in cpxControl together with any settings given
%              in cpxControl. The settings given as input are considered as
%              constants during the tuning process.
%
%   presolve   Flag controlling the CPLEX Presolve feature. The following
%              values are recognized:
%
%              0 - No special treatment of presolve (default).
%              1 - Invoke CPLEX Presolve on problem before optimization begins and
%                  create information about the changes made.
%              2 - As 1, but also returns the presolved problem, including linear
%                  constraints, objective and bounds.
%
% ------------------------------------------------------------------------------
%
% OUTPUT:
%
% x         Solution vector with decision variable values (n x 1 vector). If the
%           solution pool feature was enabled for a MIP/MIQP/MIQCP problem, x
%           may have more than one column. The first column is always the
%           incumbent IP solution.
%
% slack     Slack variables (m x 1 vector). If solution pool is active, slack
%           will have more than one column. If quadratic constraints are
%           present, the corresponding slacks are at the end of each column.
%
% v         Lagrangian multipliers (dual solution vector) (m x 1 vector)
%
% rc        Reduced costs. Lagrangian multipliers for simple bounds on x.
%
% f_k       Objective function value at optimum. If the solution pool feature
%           is enabled, f_k will be a row vector with objective values for each
%           solution. Each element in f_k thus corresponds to a column in x. 
%
% ninf      Number of infeasibilities
% sinf      Sum of infeasibilities
%
% Inform    Result of CPLEX run. See help cplexStatus.
%
% basis     Vector containing basis at solution
%
%           The first m elements contain row basis information, with the
%           following possible values for non-ranged rows:
%
%         0 associated slack/surplus/artificial variable nonbasic at value 0.0
%         1 associated slack/surplus/artificial variable basic
%
%           and for ranged rows (both upper and lower bounded)
%
%         0 associated slack/surplus/artificial variable nonbasic at its lower bound
%         1 associated slack/surplus/artificial variable basic
%         2 associated slack/surplus/artificial variable nonbasic at its upper bound
%
%           The last n elements, i.e. basis(m+1:m+n) contain column
%           basis information:
%
%         0 variable at lower bound
%         1 variable is basic
%         2 variable at upper bound
%         3 variable free and nonbasic
%
% lpiter    Number of LP iterations
%
% glnodes   Number of nodes visited
%
% confstat   Structure with extended conflict status information. This
%            output is a replica of the Prob.CPLEX.confgrps input argument
%            with the added fields 'status' and 'istat'. confstat(k).status
%            gives a text description of the status of conflict group k;
%            the corresponding istat field is the numeric value also
%            available in iconfstat(k).
%
% iconfstat  Conflict status information. For an infeasible problem where
%            at least one conflict group have been supplied in the confgrps
%            input argument, this output argument contains the status of
%            each conflict group, in the same order as given in the confgrps input.
%
%            The following values are possible:
%
%           -1  Excluded
%            0  Possible member
%            1  Possible LB
%            2  Possible UB
%            3  Member
%            4  Upper bound
%            5  Lower bound
%
%            If confstat is empty even though Conflict Refinement has been
%            requested, there was a problem in the refinement process.
%
% sa        Structure with information about the requested SA, if requested.
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
% cpxControl  Structure with non-default CPLEX parameters generated during
%             CPLEX Parameter Tuning.
%
% presolve   Structure with information about the changes CPLEX Presolve made to
%            the problem before solving. The amount of information depends on
%            the value of the cpxSettings.presolve flag.
%
%            For cpxSettings.presolve=0, this output is empty.
%
%            For cpxSettings.presolve=1, the presolve output contains six arrays,
%            pcstat, prstat, ocstat, orstat, status and objoffset.
%
%    status  Flag telling status of presolve results:
%
%         0  Problem was not presolved or no reductions were made
%         1  A presolved problem exists
%         2  The original problem was reduced to an empty problem
%
%            The remaining fields of presolve will contain useful
%            information only if status==1.
%
%    pcstat  Contains information about variables in the original
%            problem. For each element pcstat(i):
%
%        >=1 variable i corresponds to variable pcstat(i) in the presolved problem
%         -1 variable i is fixed to its lower bound
%         -2 variable i is fixed to its upper bound
%         -3 variable i is fixed to some other value
%         -4 variable i is aggregated out
%         -5 variable i is deleted or merged for some other reason
%
%    prstat  Contains information about constraints (rows) in the original
%            problem. For each element prstat(i):
%
%        >=1 row i corresponds to row prstat(i) in the presolved problem
%         -1 row i is redundant
%         -2 row i is used for aggregation
%         -3 row i is deleted for some other reason
%
%    ocstat  Contains information about variables in the presolved problem:
%            For each element ocstat(i):
%
%        >=1 variable i in the presolved problem corresponds to variable ocstat(i)
%            in the original problem.
%         -1 variable i corresponds to a linear combination of some variables in
%            the original problem
%
%    orstat  Contains information about constraints (rows) in the presolved
%            problem. For each element orstat(i):
%
%        >=1 row i in the presolved problem corresponds to row orstat(i)
%            in the original problem
%         -1 row i is created by, for example, merging two rows in the
%         original problem
%
% -----------------------------------------------------------------------
%
% NOTE! After each call of cplex, output information is defined in a
% vector called cpxRetVec, to access this global vector do
%        global cpxRetVec
%
% See the TOMLAB /CPLEX User's Guide for a description of each element

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2011 by Tomlab Optimization Inc., $Release: 12.2.0$
% Written Aug 5, 2002.    Last modified May 5, 2011.

function [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, ...
   glnodes, confstat, iconfstat, sa, cpxControl, presolve, kappa] = ...
   cplex(c, A, x_L, x_U, b_L, b_U, ...
   cpxControl, callback, PriLev, Prob, IntVars, PI, SC, SI, ...
   sos1, sos2, F, logfile, savefile, savemode, qc, ...
   confgrps, conflictFile, saRequest, basis, xIP, logcon, branchprio, ...
   branchdir, cpxSettings)

%#function cpx2cbinfo cpx2retvec
if nargin < 30
 cpxSettings = [];
if nargin <  29
 branchdir = [];
if nargin < 28
 branchprio = [];
if nargin < 27
 logcon = [];
if nargin < 26
 xIP = [];
if nargin < 25
  basis = [];
if nargin < 24
 saRequest = [];
 if nargin < 23
  conflictFile = '';
  if nargin < 22
   confgrps = [];
   if nargin < 21    
    qc = [];
    if nargin < 20
     savemode = [];
     if nargin < 19
      savefile = '';
      if nargin < 18
       logfile = '';
       if nargin < 17
        F = [];
        if nargin < 16
         sos2 = [];
         if nargin < 15
          sos1 = [];
          if nargin < 14
           SI = [];
           if nargin < 13
            SC = [];
            if nargin < 12
             PI = [];
             if nargin < 11
              IntVars = [];
              if nargin < 10
               Prob = [];
               if nargin < 9
                PriLev = [];
                if nargin < 8
                 callback = [];
                 if nargin < 7
                  cpxControl = [];
                  if nargin < 6
                   b_U = [];
                   if nargin < 5
                      error('cplex needs at least 5 arguments');
end, end, end, end, end, end, end, end, end, end, end, end, end, end , end, end
end, end, end, end, end, end, end, end, end, end

global cpxParameters cpxProblemAttributes

% Set empty to avoid old values from last run to still be used
cpxParameters        = [];
cpxProblemAttributes = [];

DEBUG=0;

n = max( [length(c),size(F,1),size(A,2),length(x_L),length(x_U)] );
if n==0, error('cplex: cannot determine problem dimension'); end

b_L = double(b_L(:));
b_U = double(b_U(:));
x_L = double(x_L(:));
x_U = double(x_U(:));

if isempty(x_L),    x_L = zeros(n,1); end
if isempty(b_U),    b_U = b_L; end
if isempty(PriLev), PriLev = 0; end
if isempty(callback) 
   callback=zeros(15,1); 
end
if length(callback)<15
   callback(15)=0;
end

[m,nA] = size(A);

% Error checking 

if nA ~= n & m > 0
   fprintf('n = %d. Number of columns in A = %d\n',n,nA);
   error('cplex: Illegal length of A');
end

if length(b_L)~=m
   fprintf('m = %d. Length of b_L = %d\n',m,length(b_L));
   error('cplex: Illegal length of b_L');
end
if length(b_U)~=m
   fprintf('m = %d. Length of b_U = %d\n',m,length(b_U));
   error('cplex: Illegal length of b_U');
end
if length(x_L)~=n
   fprintf('n = %d. Length of x_L = %d\n',n,length(x_L));
   error('cplex: Illegal length of x_L');
end
if length(x_U)~=n
   fprintf('n = %d. Length of x_U = %d\n',n,length(x_U));
   error('cplex: Illegal length of x_U');
end

if ~isempty(IntVars) | ~isempty(sos1) | ~isempty(sos2) | ~isempty(PI) | ...
   ~isempty(SC) | ~isempty(SI)
   MIP=1;
else
   MIP=0;
end

% b = ~isfinite(sum(sum(a)));

% Check f,F,A for Inf or NaN elements
if ~isfinite(sum(c))
   error('cplex: Input c contains at least one NaN or Inf element');
end
if ~isfinite(sum(sum(F)))
   error('cplex: Input F contains at least one NaN or Inf element');
end
if ~isfinite(sum(sum(A)))
   error('cplex: Input A contains at least one NaN or Inf element');
end

% Check bounds for NaN elements (Inf is OK) 
if any(isnan(x_L))
   error('cplex: Input x_L contains at least one NaN element');
end
if any(isnan(x_U))
   error('cplex: Input x_U contains at least one NaN element');
end
if any(isnan(b_L))
   error('cplex: Input b_L contains at least one NaN element');
end
if any(isnan(b_U))
   error('cplex: Input b_U contains at least one NaN element');
end

% if ~isempty(qc) & MIP==1
%   error('cplex: MILP/MIQP problems cannot have quadratic constraints');
% end

% Vector for integers indicating type
% 0    Continuous
% 1    General integer
% 2    Binary
% 3    Semi-continuous variables, 0 and real interval
% 4    Semi-integer variables, 0 and integer interval
% 5    Partially integer variables (not in CPLEX, still used now)

% Logical vector for integers
IV = zeros(n,1);

if MIP
   callback(8)=0;  % Avoid simplex log callback for MIP

   if isempty(IntVars)
      % No binary variables B or integer variables of type I
   elseif any(IntVars==0) | all(IntVars==1)
      % Assume binary logical vector given
      IV(1:length(IntVars)) = logical(IntVars);
   else
      if any(IntVars < 1 | IntVars > n)
         error('cplex: Illegal IntVars vector');
      end
      IV(IntVars)=1;
   end

   % Semi-continuous variables, o and real interval
   if ~isempty(SC)
      if any(SC < 1 | SC > n)
         error('cplex: Illegal index in SC vector');
      end
      IV(SC) = 3;
   end
   % Semi-integer variables, 0 and integer interval
   if ~isempty(SI)
      if any(SI < 1 | SI > n)
         error('cplex: Illegal index in SI vector');
      end
      IV(SI) = 4;
   end
   % Partially integer variables
   if ~isempty(PI)
      if isfield(PI,'var')
         IV(PI.var) = 5;
      else
         PI.var = [];
      end
   else
      PI.var = [];
   end

   glcolidx = find(IV);

   ngl=length(glcolidx);
   if ngl > 0
      % Identify binary variables among the general Integer variables
      ix = find(x_L(glcolidx)==0 & (x_U(glcolidx)==1) & (IV(glcolidx)==1));
      if ~isempty(ix)
         IV(glcolidx(ix)) = 2;
      end
   end

   % xIP
   if ~isempty(xIP)
      % Extend missing values at end with NaN
      if length(xIP) < n
         xIP(end+1:n) = NaN;
      end
      if ~isfield(cpxControl,'MIPSTART')
         cpxControl.MIPSTART = 1;
      end
   end
   
   if isempty(sos1) & isempty(sos2)
      settype   = [];
      setbeg    = [];
      setcolidx = [];
      setref    = [];
   else
      ns1       = length(sos1);
      ns2       = length(sos2);
      settype   = [ones(ns1,1);2*ones(ns2,1)];
      setbeg    = 1;
      setcolidx = [];
      setref    = [];
      for i=1:ns1
          if ~isfield(sos1(i),'var')
             fprintf('sos1 set %d. ',i);
             error('sos1 field var is missing');
          end
          if ~isfield(sos1(i),'row')
             fprintf('sos1 set %d. ',i);
             error('sos1 field row is missing');
          end
          ix        = sos1(i).var;
          if ~(all(ix >= 1 & ix <= n))
             fprintf('sos1 set %d. ',i);
             fprintf('\n');
             error('cplex: Illegal sos1 input variable vector');
          end
          row       = sos1(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos1 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('cplex: Illegal sos1 row data');
          end
          k         = length(ix);
          setbeg    = [setbeg; setbeg(length(setbeg))+k];
          setcolidx = [setcolidx; ix(:)];
          if row==0
             key = full(c(ix));
             %key = full(c(ix)) + 0.05*rand(length(ix),1);
          else
             key = full(A(row,ix))';
             %key = full(A(row,ix))' + 0.05*rand(length(ix),1);
          end
          % Sort the ordering key, to get unique sequence of variables
          [zz,skey]=sort(key);
          setref = [setref; skey(:)];
      end
      if DEBUG
         xprinti(setbeg,'setbeg');
         xprinti(setref,'setref');
         xprinti(setcolidx,'setcolidx');
         pause
      end
      for i=1:ns2
          if ~isfield(sos2(i),'var')
             fprintf('sos2 set %d. ',i);
             error('sos2 field var is missing');
          end
          if ~isfield(sos2(i),'row')
             fprintf('sos2 set %d. ',i);
             error('sos2 field row is missing');
          end
          ix        = sos2(i).var;
          if ~(all(ix >= 1 & ix <= n))
             fprintf('sos2 set %d. ',i);
             fprintf('\n');
             error('cplex: Illegal sos2 input variable vector');
          end
          row       = sos2(i).row;
          if ~(row >= 0 & row <= m)
             fprintf('sos2 set %d. ',i);
             fprintf('Illegal row number  %d.',row);
             fprintf('\n');
             error('cplex: Illegal sos2 row data');
          end
          k         = length(ix);
          setbeg    = [setbeg; setbeg(length(setbeg))+k];
          setcolidx = [setcolidx; ix(:)];
          if row==0
             key = full(c(ix));
          else
             key = full(A(row,ix))';
          end
          % Sort the ordering key, to get unique sequence of variables
          [zz,skey]=sort(key);
          setref = [setref; skey(:)];
      end
   end
else
   gltype    = []; glcolidx  = []; gllim     = []; settype   = [];
   setbeg    = []; setcolidx = []; setref    = [];
end

% Check that Prob.P is set.
if isempty(Prob)
    Prob = struct('P',1);
end
if ~isfield(Prob,'P')
    Prob.P=1;
end

if ~isempty(branchprio)
    if length(branchprio) ~= n | any(branchprio < 0)
        error('branchprio incorrectly defined. Check size and that entries are >= 0')
    end
end

if ~isempty(branchdir)
    if length(branchdir) ~= n
        error('branchdir incorrectly defined. Check size')
    end
end

if any(callback)
    % Define fields in Prob for use in callback m-files
    Prob.A    = sparse(A);
    Prob.QP.c = c;
    Prob.QP.F = sparse(F);
    Prob.QP.qc = qc;
    Prob.x_L  = x_L;
    Prob.x_U  = x_U;
    Prob.b_L  = b_L;
    Prob.b_U  = b_U;
    Prob.iv        = IV;
    Prob.qstype    = settype;
    Prob.msstart   = setbeg;
    Prob.mscols    = setcolidx;
    Prob.dref      = setref;
end

% The cplex MEX has 1E12 as default, but is changed if Prob.BIG is noempty
BIG=DefPar(Prob,'BIG',1E12);

try
   [x, slack, v, rc, f_k, ninf, sinf, Inform, basis, lpiter, glnodes,...
      iconfstat, sa, cpxControl, presolve, retvec ] ...
      = cplexmex(double(c), sparse(F), sparse(A), callback , x_L, x_U, b_L, b_U,...
      BIG, cpxSettings, PriLev, Prob, ...
      IV, settype, setbeg, setcolidx, setref, cpxControl, branchprio, branchdir, ...
      logfile, savefile, savemode, qc, confgrps, conflictFile, saRequest, ...
      basis, xIP, logcon );
catch
   l=lasterror;
   if(strcmp(l.identifier,'MATLAB:invalidMEXFile'))
      tomlabsharederror;
   else
      rethrow(l);
   end
end

% Condition number statistics
try
   kappa = struct('Attention',retvec(82),'IllPosed',retvec(83),...
    'Max',retvec(84),'Stable',retvec(85),...
    'Suspicious', retvec(86),'Unstable',retvec(87));
catch
   kappa = struct('Attention',NaN, 'IllPosed',NaN,...
    'Max', NaN, 'Stable', NaN,...
    'Suspicious', NaN, 'Unstable', NaN );
end

if ~isempty(iconfstat)
   confstat = confgrps;
   for k=1:length(iconfstat)
      confstat(k).istat  = iconfstat(k);
      switch(iconfstat(k))
         case -1, s = 'excluded';  
         case  0, s = 'possible'; 
         case  1, s = 'possible LB'; 
         case  2, s = 'possible UB'; 
         case  3, s = 'member'; 
         case  4, s = 'LB';
         case  5, s = 'UB';
         otherwise, s = 'unknown';
      end
      confstat(k).status = s;
   end
else
   confstat = [];
end

% MODIFICATION LOG:
%
% 020805 fhe  Modified from xpress.m  to be used with CPLEX
% 020806 hkh  Revision for CPLEX use, change integer variable handling
% 020922 hkh  Revision of comments
% 030731 ango Correction of comments on Inform
% 030925 ango New dump mode (LP) added; logfile feature added
% 030930 ango Removed previous handling of problem dump
% 031125 ango CPLEX 9.0 adaptations. Error codes reviewed
% 040101 hkh  Wrong computation of n, must use x_L,x_U, not b_L and b_U
% 040803 med  Added pragmas for MATLAB Compiler
% 040818 fhe  Added iis and sa and some help text for iis and sa.
% 040825 ango Slight change for sa structure fields
% 041213 hkh  Set BIG hard coded as 1E12, if not Prob.BIG is set
% 041213 hkh  Add comments about BIG and inf
% 050128 med  Function DefPar added
% 050201 ango Basis support added
% 050202 med  Help updated
% 050414 ango mipstart added, 9.1 ready
% 050504 ango Change mipstart, use xIP instead
% 050707 med  try-catch statement added to find install problem
% 050922 med  Updated error due to install problem
% 060129 ango CPLEX 10.0 adaptations
% 060207 ango try-catch updated
% 060818 med  Force c, x_L, x_U, b_L and b_U to doubles
% 061207 ango Check for inf/nan in A, F, c and bounds vectors
% 070222 hkh  IntVars format revised
% 070526 med  mlint check
% 070601 ango Incumbent callback added
% 070611 med  CPLEX return codes corrected
% 070802 ango User cut callback added
% 071116 med  callback print out removed
% 071116 med  branchprio and branchdir added
% 071213 ango Parameter Tuning added
% 080204 ango Presolve feature added
% 080207 med  Help corrected
% 080220 ango Change obsolete parameter MIPSTART into ADVIND
% 090804 med  Removed version check
