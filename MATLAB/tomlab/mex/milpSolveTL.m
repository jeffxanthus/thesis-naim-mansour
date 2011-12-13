% TOMLAB MILPSOLVE 5.1 (MI)LP Solver
%
% milpSolve is a mixed-integer linear programming solver based on the
% open source LP_SOLVE licensed under GNU LGPL. Complete source code for the
% solver engine and associated BFP packages is available here:
%
% http://groups.yahoo.com/group/lp_solve
%
% -------------------------------------------------------------------------
%
% function Result = milpSolveTL(Prob) called by one of the following options:
%
% Result = tomRun('MILPSOLVE', Prob, 1); or
% Prob = ProbCheck(Prob, 'MILPSOLVE');
% Result = milpSolveTL(Prob);
% PrintResult(Result,1);
%
% INPUT:
%
%   Prob        Problem structure in TOMLAB format. Use lpAssign or
%               mipAssign to define the Prob structure.
%
%   Fields used in input structure Prob:
%
%   x_L, x_U    Lower and upper bounds on variables. (Must be dense)
%   b_L, b_U    Lower and upper bounds on linear constraints. (Must be dense)
%   A           Linear constraint matrix. (Sparse or dense)
%
%   QP.c        Linear objective function coefficients, size n x 1
%
%   BIG         Definition of infinity. Default is 1e30
%
%   LargeScale  Defines if milpSolveTL will convert the A matrix to a sparse
%               matrix or not.
%               LargeScale != 0 - sparse
%               LargeScale  = 0 - dense
%               Default is to use the A matrix just as it is defined
%
%   PriLevOpt   Specifies the printlevel that will be used by MILPSOLVE.
%
%           0 (NONE)      No outputs
%           1 (NEUTRAL)   Only some specific debug messages in debug print
%           routines are reported.
%           2 (CRITICAL)  Only critical messages are reported. Hard errors
%           like instability, out of memory.
%           3 (SEVERE)    Only severe messages are reported. Errors.
%           4 (IMPORTANT) Only important messages are reported. Warnings
%           and Errors.
%           5 (NORMAL)    Normal messages are reported.
%           6 (DETAILED)  Detailed messages are reported. Like model size,
%           continuing B&B improvements.
%           7 (FULL)      All messages are reported. Useful for debugging
%           purposes and small models.
%
%           Default print level is 0, no outputs.
%           PriLevOpt < 0 is interpreted as 0, and larger than 7 is interpreted
%           as 7.
%
%   MaxCPU      Maximal CPU Time (in seconds) to be used by MILPSOLVE,
%               stops with best point found.
%
% Fields used in Prob.MILPSOLVE (Structure with MILPSOLVE specific parameters)
%
%   ANTI_DEGEN  Binary vector. If empty, no anti-degeneray handling is
%               applied. If the length (i) of the vector is less than 8 elements,
%               only the i first modes are considered. Also if i is longer than
%               8 elements, the elements after element 8 are ignored.
%
%           ANTI_DEGEN specifies if special handling must be done to reduce
%           degeneracy/cycling while solving. Setting this flag can avoid
%           cycling, but can also increase numerical instability.
%
%           ANTIDEGEN_FIXEDVARS   != 0  Check if there are equality slacks in
%           the basis and try to drive them out in order to reduce chance of
%           degeneracy in Phase 1
%           ANTIDEGEN_COLUMNCHECK != 0
%           ANTIDEGEN_STALLING    != 0
%           ANTIDEGEN_NUMFAILURE  != 0
%           ANTIDEGEN_LOSTFEAS    != 0
%           ANTIDEGEN_INFEASIBLE  != 0
%           ANTIDEGEN_DYNAMIC     != 0
%           ANTIDEGEN_DURINGBB    != 0
%
%   basis   If empty or erroneous, default basis is used. Default start
%           base is the all slack basis (the default simplex starting basis).
%
%           Prob.MILPSOLVE.basis stores the basic variables. If an element is
%           less then zero then it means on lower bound, else on upper bound.
%           Element 0 of the array is unused. The default initial basis is
%           bascolumn[x] = -x. By MILPSOLVE convention, a basic variable is
%           always on its lower bound, meaning that basic variables is
%           always represented with a minus sign.
%
%           When a restart is done, the basis vector must be assigned a correct
%           starting basis and Prob.WarmStart must be set to 1.
%
%   BASIS_CRASH Specifies which basis crash mode will used.
%
%           When no base crash is done (the default), the initial basis
%           from which MILPSOLVE starts to solve the model is the basis
%           containing all slack or artificial variables that is
%           automatically associates with each constraint.
%
%           When base crash is enabled, a heuristic ``crash procedure''
%           is executed before the first simplex iteration to quickly
%           choose a basis matrix that has fewer artificial variables.
%           This procedure tends to reduce the number of iterations to
%           optimality since a number of iterations are skipped.
%           MILPSOLVE starts iterating from this basis until optimality.
%
%              BASIS_CRASH != 2  - No basis crash
%              BASIS_CRASH  = 2  - Most feasible basis
%
%           Default is no basis crash.
%
%   BB_DEPTH_LIMIT Sets the maximum branch-and-bound depth. This value makes
%           sense only if there are integer, semi-continuous or SOS
%           variables in the model so that the branch-and-bound algorithm
%           is used to solve the model. The branch-and-bound algorithm will not
%           go deeper than this level. When BB_DEPTH_LIMIT i set to 0 then there
%           is no limit to the depth. The default value is -50. A positive value
%           means that the depth is absolute. A negative value means a relative
%           B&B depth. The "order" of a MIP problem is defined to be 2
%           times the number of binary variables plus the number of SC
%           and SOS variables. A relative value of -x results in a maximum
%           depth of x times the order of the MIP problem.
%
%   BB_FLOOR_FIRST Specifies which branch to take first in branch-and-bound
%           algorithm. Default value is 1.
%
%           BB_FLOOR_FIRST = 0 (BRANCH_CEILING)   Take ceiling branch first
%           BB_FLOOR_FIRST = 1 (BRANCH_FLOOR)     Take floor branch first
%           BB_FLOOR_FIRST = 2 (BRANCH_AUTOMATIC) MILPSOLVE decides which branch
%           being taken first
%
%   BB_RULE     Specifies the branch-and-bound rule. Default value is 0.
%
%           BB_RULE = 0 (NODE_FIRSTSELECT)        Select lowest indexed
%           non-integer column
%           BB_RULE = 1 (NODE_GAPSELECT)          Selection based on distance
%           from the current bounds
%           BB_RULE = 2 (NODE_RANGESELECT)        Selection based on the
%           largest current bound
%           BB_RULE = 3 (NODE_FRACTIONSELECT)     Selection based on largest
%           fractional value
%           BB_RULE = 4 (NODE_PSEUDOCOSTSELECT4)  Simple, unweighted pseudo
%           -cost of a variable
%           BB_RULE = 5 (NODE_PSEUDONONINTSELECT) This is an extended pseudo
%           -costing strategy based on minimizing the number of integer infeasibilities.
%           BB_RULE = 6 (NODE_PSEUDORATIOSELECT)  This is an extended pseudo
%           -costing strategy based on maximizing the normal pseudo-cost
%           divided by the number of infeasibilities. Effectively, it is
%           similar to (the reciprocal of) a cost/benefit ratio.
%           BB_RULE = 7 (NODE_USERSELECT)
%
%   BB_RULE_ADD Additional values for the BB_RULE. BB_RULE is a vector.
%           If the length i of the vector is less than 10 elements, only
%           the i first modes are considered. Also if i is longer than
%           10 elements, the elements after element 10 is ignored.
%
%           BB_RULE_ADD(1) != 0  (NODE_WEIGHTREVERSEMODE)
%           BB_RULE_ADD(2) != 0  (NODE_BRANCHREVERSEMODE) In case when
%           BB_FLOOR_FIRST is BRANCH_AUTOMATIC, select the opposite
%           direction (lower/upper branch) that BRANCH_AUTOMATIC had chosen.
%           BB_RULE_ADD(3) != 0  (NODE_GREEDYMODE)
%           BB_RULE_ADD(4) != 0  (NODE_PSEUDOCOSTMODE)
%           BB_RULE_ADD(5) != 0  (NODE_DEPTHFIRSTMODE)    Select the node
%           that has already been selected before the number of times
%           BB_RULE_ADD(6) != 0  (NODE_RANDOMIZEMODE)
%           BB_RULE_ADD(7) != 0  (NODE_DYNAMICMODE)       When NODE_DEPTHFIRSTMODE
%           is selected, switch off this mode when a first solution is found.
%           BB_RULE_ADD(8) != 0  (NODE_RESTARTMODE)
%           BB_RULE_ADD(9) != 0  (NODE_BREADTHFIRSTMODE)  Select the node
%           that has been selected before the fewest number of times or not at all
%           BB_RULE_ADD(10) != 0 (NODE_AUTOORDER)
%
%   BFP     Defines which Basis Factorization Package that will be used
%           by MILPSOLVE.
%
%           BFP = 0 : LUSOL
%           BFP = 1 : built in etaPHI from LP_SOLVE v3.2
%           BFP = 2 : Additional etaPHI
%           BFP = 3 : GLPK
%
%           Default BFP is LUSOL.
%
%   BREAK_AT_FIRST Specifies if the branch-and-bound algorithm stops at the first
%           found solution (BREAK_AT_FIRST != 0) or not (BREAK_AT_FIRST = 0).
%           Default is not to stop at the first found solution.
%
%   BREAK_AT_VALUE Specifies if the branch-and-bound algorithm stops when the
%           object value is better than a given value. The default value
%           is (-) infinity.
%
%   EPAGAP      Specifies the absolute MIP gap tolerance for the branch and bound
%               algorithm. This tolerance is the difference between the best-found
%               solution yet and the current solution. If the difference is smaller
%               than this tolerance then the solution (and all the sub-solutions)
%               is rejected. The default value is 1e-9.
%
%   EPGAP       Specifies the relative MIP gap tolerance for the branch and
%               bound algorithm. The default value is 1e-9.
%
%   EPSB        Specifies the value that is used as a tolerance for the Right
%               Hand Side (RHS) to determine whether a value should be considered
%               as 0. The default epsb value is 1.0e-10
%
%   EPSD        Specifies the value that is used as a tolerance for reduced
%               costs to determine whether a value should be considered as 0.
%               The default epsd value is 1e-9. If EPSD is empty, EPSD is read from
%               Prob.optParam.eps_f.
%
%   EPSEL       Specifies the value that is used as a tolerance for rounding
%               values to zero. The default epsel value is 1e-12.
%
%   EPSINT      Specifies the tolerance that is used to determine whether a
%               floating-point number is in fact an integer. The default value
%               for epsint is 1e-7. Changing this tolerance value can result
%               in faster solving times, but the solution is less integer.
%
%   EPSPERTURB  Specifies the value that is used as perturbation scalar for
%               degenerative problems. The default epsperturb value is 1e-5
%
%   EPSPIVOT    Specifies the value that is used as a tolerance pivot element
%               to determine whether a value should be considered as 0.
%               The default epspivot value is 2e-7
%
%   IMPROVEMENT_LEVEL Specifies the iterative improvement level.
%
%              IMPROVEMENT_LEVEL = 0  (IMPROVE_NONE)    improve none
%              IMPROVEMENT_LEVEL = 1  (IMPROVE_FTRAN)   improve FTRAN
%              IMPROVEMENT_LEVEL = 2  (IMPROVE_BTRAN)   improve BTRAN
%              IMPROVEMENT_LEVEL = 3  (IMPROVE_SOLVE)   improve FTRAN + BTRAN.
%              IMPROVEMENT_LEVEL = 4  (IMPROVE_INVERSE) triggers automatic
%           inverse accuracy control in the dual simplex, and when an error
%           gap is exceeded the basis is reinverted
%
%               Choice 1,2,3 should not be used with LP_SOLVE 5.1.1.3, because
%               of problems with the solver. Default is 0.
%
%   LoadFile    File that contains the model. If LoadFile is a nonempty string
%               which corresponds to actual file, then the model is read from
%               this file rather than from the Prob struct.
%
%   LoadMode    1 - LP   - MILPSOLVE LP format
%               2 - MPS  - MPS format
%               3 - FMPS - Free MPS format
%
%               A default value for this field does not exist. Both LoadFile
%               and LoadMode must be set if a problem will be loaded.
%
%               If there is something wrong with LoadMode or LoadFile, an error
%               message will be printed and MILPSOLVE will be terminated. Leave
%               LoadMode and LoadFile empty if the problem not will be loaded
%               from file.
%
%   LogFile     Name of file to print MILPSOLVE log on.
%
%
%   MAXIMIZE    If MAXIMIZE != 0, MILPSOLVE is set to maximize the objective
%               function, default is to minimize.
%
%   MAX_PIVOT   Sets the maximum number of pivots between a re-inversion of
%               the matrix. Default is 42.
%
%   NEG_RANGE   Specifies the negative value below which variables are split
%               into a negative and a positive part. This value must always be
%               zero or negative. If a positive value is specified, then 0 is taken.
%               The default value is -1e6.
%
%   PRESOLVE    Vector containing possible presolve options. If the length i
%               of the vector is less than 7 elements, only the i first modes
%               are considered. Also if i is longer than 7 elements, the
%               elements after element 7 is ignored.
%
%               PRESOLVE(1) != 0 (PRESOLVE_ROWS)      Presolve rows
%               PRESOLVE(2) != 0 (PRESOLVE_COLS)      Presolve columns
%               PRESOLVE(3) != 0 (PRESOLVE_LINDEP)    Eliminate linearly dependent rows
%               PRESOLVE(4) != 0 (PRESOLVE_SOS)       Convert constraints to SOSes
%               (only SOS1 handled)
%               PRESOLVE(5) != 0 (PRESOLVE_REDUCEMIP) If the phase 1 solution
%               process finds that a constraint is redundant then this constraint
%               is deleted.
%               PRESOLVE(6) != 0 (PRESOLVE_DUALS)     Calculate duals
%               PRESOLVE(7) != 0 (PRESOLVE_SENSDUALS) Calculate sensitivity if
%               there are integer variables
%
%               Default is not to do any presolve.
%
%   PRICING_RULE The pricing rule can be one of the following rules.
%
%            PRICING_RULE = 0  Select first (PRICER_FIRSTINDEX)
%            PRICING_RULE = 1  Select according to Dantzig (PRICER_DANTZIG)
%            PRICING_RULE = 2  Devex pricing from Paula Harris (PRICER_DEVEX)
%            PRICING_RULE = 3  Steepest Edge (PRICER_STEEPESTEDGE)
%
%   PRICING_MODE Additional pricing settings, any combination of the
%             modes below. This is a binary vector. If the length i of the
%             vector is less than 7 elements, only the i first modes are
%             considered. Also if i is longer than 7 elements, the elements
%             after element 7 is ignored.
%
%            PRICE_PRIMALFALLBACK  != 0 In case of Steepest Edge, fall back to
%            DEVEX in primal.
%            PRICE_MULTIPLE        != 0 Preliminary implementation of the
%            multiple pricing scheme. This means that attractive candidate
%            entering columns from one iteration may be used in the subsequent
%            iteration, avoiding full updating of reduced costs. In the current
%            implementation, MILPSOLVE only reuses the 2nd best entering column
%            alternative.
%            PRICE_PARTIAL         != 0 Enable partial pricing
%            PRICE_ADAPTIVE        != 0 Temporarily use First Index if cycling
%            is detected
%            PRICE_RANDOMIZE       != 0 Adds a small randomization effect to
%            the selected pricer
%            PRICE_LOOPLEFT        != 0 Scan entering/leaving columns left
%            rather than right
%            PRICE_LOOPALTERNATE   != 0 Scan entering/leaving columns
%            alternatingly left/right
%
%            Default basis is PRICE_DEVEX combined with PRICE_ADAPTIVE.
%
%   sa       Struct containing information of the sensitivity analysis (SA)
%            MILPSOLVE will perform.
%
%            sa.obj =! 0 Perform sensitivity analysis on the objective function
%            sa.obj =  0 Do not perform sensitivity analysis on the objective function
%            sa.rhs =! 0 Perform sensitivity analysis on the right hand sides.
%            sa.rhs =  0 Do not perform sensitivity analysis on the right hand sides.
%
%   SaveFileAfter Name of a file to save the MILPSOLVE object after presolve.
%            The name must be of type string (char),
%            Example: Prob.MILPSOLVE.SaveFileAfter = 'save2'
%            If the type is not char SaveFileBefore is set to save2.[file_extension].
%
%   SaveFileBefore Name of a file to save the MILPSOLVE object before presolve.
%            The name must be of type string (char),
%            Example: Prob.MILPSOLVE.SaveFileBefore = 'save1'
%            If the type is not char SaveFileBefore is set to save1.[file_extension].
%
%   SaveMode 1 - LP   - MILPSOLVE LP format
%            2 - MPS  - MPS format
%            3 - FMPS - Free MPS format
%            If empty, the default format LP is used.
%
%   SCALE_LIMIT Sets the relative scaling convergence criterion to the absolute
%             value of SCALE_LIMIT for the active scaling mode.
%             The integer part of SCALE_LIMIT specifies the maximum number
%             of iterations. Default is 5.
%
%   SCALING_ALG Specifies which scaling algorithm will be used by MILPSOLVE.
%
%            0 No scaling algorithm
%            1 (SCALE_EXTREME) Scale to convergence using largest
%            absolute value
%            2 (SCALE_RANGE) Scale based on the simple numerical range
%            3 (SCALE_MEAN) Numerical range-based scaling
%            4 (SCALE_GEOMETRIC) Geometric scaling
%            7 (SCALE_CURTISREID) Curtis-reid scaling
%
%            Default is 0, no scaling algorithm.
%
%   SCALING_ADD Vector containing possible additional scaling parameters.
%            If the length (i) of the vector is less than 7 elements, only
%            the i first modes are considered. Also if i is longer than 7
%            elements, the elements after element 7 is ignored.
%
%            SCALING_ADD != 0  (SCALE_QUADRATIC)
%            SCALING_ADD != 0  (SCALE_LOGARITHMIC)  Scale to convergence using
%            logarithmic mean of all values
%            SCALING_ADD != 0  (SCALE_USERWEIGHT)   User can specify scalars
%            SCALING_ADD != 0  (SCALE_POWER2)       also do Power scaling
%            SCALING_ADD != 0  (SCALE_EQUILIBRATE)  Make sure that no scaled
%            number is above 1
%            SCALING_ADD != 0  (SCALE_INTEGERS)     Also scaling integer variables
%            SCALING_ADD != 0  (SCALE_DYNUPDATE)    Dynamic update
%
%            Default is 0, no additional mode.
%
%            Settings SCALE_DYNUPDATE is a way to make sure that scaling
%            factors are recomputed. In that case, the scaling factors are
%            recomputed also when a restart is done.
%
%   SIMPLEX_TYPE Sets the desired combination of primal and dual simplex algorithms.
%
%            5  (SIMPLEX_PRIMAL_PRIMAL) Phase1 Primal, Phase2 Primal
%            6  (SIMPLEX_DUAL_PRIMAL)   Phase1 Dual, Phase2 Primal
%            9  (SIMPLEX_PRIMAL_DUAL)   Phase1 Primal, Phase2 Dual
%            10 (SIMPLEX_DUAL_DUAL)      Phase1 Dual, Phase2 Dual
%
%            Default is SIMPLEX_DUAL_PRIMAL (6).
%
%   SOLUTION_LIMIT Sets the solution number that will be returned. This value
%            is only considered if there are integer, semi-continuous or
%            SOS variables in the model so that the branch-and-bound
%            algorithm is used. If there are more solutions with the same
%            objective value, then this number specifies which solution
%            must be returned. Default is 1.
%
%   sos       List of structs containing data about Special Ordered Sets (SOS)
%             See below for further description.
%
%
% Fields used in Prob.MIP:
%
%   IntVars   Defines which variables are integers, of general type I or
%             binary type B Variable indices should be in the range [1,...,n].
%
%             IntVars is a logical vector ==> x(find(IntVars > 0)) are integers
%
%             IntVars is a vector of indices ==> x(IntVars) are integers
%             (if [], then no integers of type I or B are defined)
%             variables with x_L=0 and x_U=1, is are set to binary.
%             It is possible to combine integer and semi-continuous type to
%             obtain the semi-integer type.
%
%   fIP       This parameter specifies the initial "at least
%             better than" guess for objective function. This is only used
%             in the branch-and-bound algorithm when integer variables exist
%             in the model. All solutions with a worse objective value than
%             this value are immediately rejected.
%             The default is infinity.
%
%   SC        A vector with indices for variables of type semi-continuous
%             (SC), a logical vector or a scalar (see MIP.IntVars).
%             A semi-continuous variable i takes either the value 0 or some
%             value in the range [x_L(i), x_U(i)]. It is possible to combine
%             integer and semi-continuous type to obtain the semi-integer type.
%
%   sos1      List of structures defining the Special Ordered Sets of Order One (SOS1).
%             For SOS1 set k,
%             sos1(k).var is a vector of indices for variables of type SOS1 in set k,
%             sos1(k).row is the priority of SOS k in the set of SOS1 and
%             sos1(k).weight is a vector of the same length as sos1(k).var and it
%             describes the order MILPSOLVE will weight the variables in SOS1 set k.
%
%             a low number of a row and a weight means high priority.
%
%   sos2      List of n structures defining the Special Ordered Sets (SOS)
%             of Order Two (SOS2). (see MIP.sos1)
%
% OUTPUT:
%
%   Result      Structure containing information about solved problem
%
% Fields used in output structure Result:
%
%   x_k       Optimal solution (or some other solution if optimum could not
%             been found)
%
%   f_k       Optimal objective value.
%
%   v_k       [rc; duals]. If Reduced cost and dual variables are not available,
%             then v_k is empty.
%
%   ExitFlag  TOMLAB information parameter.
%
%   ExitText  Status text from MILPSOLVE.
%
%   Inform    MILPSOLVE information parameter.
%
%   Iter      The total number of nodes processed in the branch-and-bound algorithm.
%             Is only applicable if the model contains integer variables.
%
%             In the case of an LP model Result.Iter contains the number of
%             iterations. This is however not documented.
%
%   MinorIter The total number of Branch-and-bound iterations.
%             When the problem is LP, MinorIter equals Result.Iter.
%
%   MILPSOLVE.basis Optimal basis, on the format described above under
%   Prob.MILPSOLVE.basis.
%
%   MILPSOLVE.MaxLevel The deepest Branch-and-bound level of the last solution.
%            Is only applicable if the model contains integer variables.
%
%   MILPSOLVE.sa Structure containing information about the sensitivity analysis.
%
%   xState     State of each variable
%              0 - free variable,
%              1 - variable on lower bound,
%              2 - variable on upper bound,
%              3 - variable is fixed, lower bound = upper bound
%
%   bState     State of each linear constraint
%              0 - Inactive constraint
%              1 - Linear constraint on lower bound
%              2 - Linear constraint on upper bound
%              3 - Linear equality constraint
%
%
% Fields used in Result.MILPSOLVE.sa
%   obj.status
%             1  successful
%             0  SA not requested
%            -1  Error: error from MILPSOLVE
%            -3  no SA available
%
%   obj.lower An array that will contain the values of the lower
%             limits on the objective function.
%
%   obj.upper An array that will contain the values of the upper
%             limits on the objective function.
%
%   rhs.status see MILPSOLVE.sa.objStatus.
%
%   rhs.lower An array that will contain the values of the lower
%             limits on the RHS.
%
%   rhs.upper An array that will contain the values of the upper
%             limits on the RHS.

% Johan Holmgren, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com
% Copyright (c) 2005-2009 by Tomlab Optimization Inc., $Release: 7.2.0$
% Written Feb 25, 2005. Last modified Aug 4, 2009.

function Result = milpSolveTL(Prob)

%#function lp_f lp_g lp_H

if nargin < 1
    error('milpSolveTL needs the Prob structure as input');
end

Prob.solvType = 7; % MILP solver
Prob = iniSolveMini(Prob);

Prob.MILPSOLVE = DefPar(Prob,'MILPSOLVE',[]);
Prob.BIG = DefPar(Prob, 'BIG', 1e30);

if isfield(Prob.MILPSOLVE, 'basis')
    if(strcmp(class(Prob.MILPSOLVE.basis),'int32') == 0)
        Prob.MILPSOLVE.basis = int32(Prob.MILPSOLVE.basis);
    end
else
    Prob.MILPSOLVE.basis = [];
end

[bl, bu, n, m] = defblbu(Prob, Prob.BIG);
nTot=n+m;
A = [Prob.QP.c(:)' ; Prob.A];

if( ~isempty(Prob.LargeScale))
    if(Prob.LargeScale == 0)
        A = full(A);
    else
        A = sparse(A);
    end
end

% The following case will give an error message in the mex:
% IntVars is given as a logical vector, but not of full length Prob.N
% But better to allow this, and avoid wasted CPU and memory here.
Prob.MIP.IntVars = DefPar(Prob.MIP,'IntVars',[]);

Prob.MIP.sos1    = DefPar(Prob.MIP,'sos1',[]);
Prob.MIP.sos2    = DefPar(Prob.MIP, 'sos2', []);
Prob.MIP.SC      = DefPar(Prob.MIP,'SC',[]);

Prob.MILPSOLVE.SaveFileBefore = DefPar(Prob.MILPSOLVE,'SaveFileBefore',[]);
Prob.MILPSOLVE.SaveFileAfter  = DefPar(Prob.MILPSOLVE,'SaveFileAfter',[]);
Prob.MILPSOLVE.SaveMode       = DefPar(Prob.MILPSOLVE,'SaveMode',[]);
Prob.MILPSOLVE.LogFile        = DefPar(Prob.MILPSOLVE,'LogFile',[]);
Prob.MILPSOLVE.LoadFile       = DefPar(Prob.MILPSOLVE,'LoadFile',[]);
Prob.MILPSOLVE.LoadMode       = DefPar(Prob.MILPSOLVE,'LoadMode',[]);
Prob.MILPSOLVE.EPSD           = DefPar(Prob.MILPSOLVE,'EPSD', []);

Result=ResultDef(Prob);
Result.Solver='milpSolve';

Prob.MILPSOLVE.SIMPLEX_TYPE = DefPar(Prob.MILPSOLVE,'SIMPLEX_TYPE',[]);
if isempty(Prob.MILPSOLVE.SIMPLEX_TYPE)
    Result.SolverAlgorithm='Phase1 Dual Simplex, Phase2 Primal Simplex';
elseif Prob.MILPSOLVE.SIMPLEX_TYPE == 5
    Result.SolverAlgorithm='Phase1 Primal Simplex, Phase2 Primal Simplex';
elseif Prob.MILPSOLVE.SIMPLEX_TYPE == 9
    Result.SolverAlgorithm='Phase1 Primal Simplex, Phase2 Dual Simplex';
elseif Prob.MILPSOLVE.SIMPLEX_TYPE == 10
    Result.SolverAlgorithm='Phase1 Dual Simplex, Phase2 Dual Simplex';
else
    Result.SolverAlgorithm='Phase1 Dual Simplex, Phase2 Primal Simplex';
end

if(isempty(Prob.MILPSOLVE.EPSD))
   Prob.MILPSOLVE.EPSD = Prob.optParam.eps_f;
end

try
[Result.x_k, Result.f_k, dual, rc, Result.MILPSOLVE.sa, Result.MILPSOLVE.basis,...
        Result.Inform, Result.ExitText, Result.MinorIter,...
        Result.Iter, Result.MILPSOLVE.MaxLevel]...
    = milpSolve(Prob, A, bl(1:n), bu(1:n), bl(n+1:nTot), bu(n+1:nTot), ...
    Prob.MIP.IntVars, Prob.MIP.SC, Prob.MIP.sos1, ...
    Prob.MIP.sos2, Prob.MILPSOLVE.basis, Prob.PriLevOpt, ...
    Prob.MILPSOLVE.LoadFile, Prob.MILPSOLVE.LoadMode, ...
    Prob.MILPSOLVE.SaveFileBefore, Prob.MILPSOLVE.SaveFileAfter, ...
    Prob.MILPSOLVE.SaveMode, Prob.MILPSOLVE.LogFile);
catch
   l=lasterror;
   if(strcmp(l.identifier,'MATLAB:invalidMEXFile'))
      tomlabsharederror;
   else
      rethrow(l);
   end
end
    
Result.v_k = [rc(:); dual(:)];

if m > 0
    Result.Ax = Prob.A * Result.x_k(:);
    Result = StateDef(Result, Result.x_k(1:n), Result.Ax, [], ...
        Prob.optParam.xTol, Prob.optParam.bTol, [], bl, bu, 1);
else
    Result = StateDef(Result, Result.x_k(1:n),[],[], Prob.optParam.xTol,[],[], ...
        bl, bu, 1);
end

switch Result.Inform
    case -2
        Result.ExitFlag = 10;
        % ExitText = 'Out of memory';
    case 0
        Result.ExitFlag = 0;
        % ExitText = 'Optimal solution found';
    case 1
        Result.ExitFlag = 1;
        % ExitText = 'Suboptimal solution';
    case 2
        Result.ExitFlag = 4;
        % ExitText = 'Infeasible model';
    case 3
        Result.ExitFlag = 2;
        %  ExitText = 'Unbounded solution';
    case 4
        Result.ExitFlag = 0;
        %  ExitText = 'Degenerate solution';
    case 5
        Result.ExitFlag = 3;
        %  ExitText = 'Numerical failure';
    case 6
        Result.ExitFlag = 1;
        %  ExitText = 'User abort';
    case 7
        Result.ExitFlag = 1;
        %  ExitText = 'Timeout';
    case 10
        Result.ExitFlag = 3;
        %  ExitText = 'Branch and bound failed';
    case 11
        Result.ExitFlag = 11;
        %  ExitText = 'Branch and bound stopped'
    case 12
        Result.ExitFlag = 0;
        %   ExitText = 'Feasible branch and bound solution';
    case 13
        Result.ExitFlag = 4;
        %  ExitText = 'No feasible branch and bound solution';
    otherwise
        Result.ExitFlag = 3;
        %  ExitText =  'Unknown return code'
end

Result.x_0=[];
Result.f_0=0;

if ~isempty(Prob.QP.c)
    Result.g_k=Prob.QP.c;
else
    Result.g_k=[];
end

Result=endSolveMini(Prob,Result);

% MODIFICATION LOG:
%
% 050324 joho Written
% 050401 joho Changes in documentation
% 050422 ango Change name to milpSolve
% 050914 med  try-catch statement added for errors
% 050922 med  Updated installation error
% 060217 med  Updated mex error catch
% 061005 med  Minor help updates
% 070222 hkh  Revise IntVars handling, use new format
% 070323 med  MaxCPU parameter added to help
% 080606 med  Switched to iniSolveMini
% 080607 hkh  Switched to endSolveMini
% 090804 med  Removed version check