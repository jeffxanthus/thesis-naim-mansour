% amplAssign implements the TOMLAB format for:
%
% LP    - Linear Programming.
% MILP  - Mixed Integer Linear Programming.
% QP    - Quadratic Programming.
% MIQP  - Mixed Integer Quadratic Programming.
% NLP   - NonLinear Programming.
% MINLP - Mixed Integer NonLinear Programming.
% LCP   - Linear Complementarity Problem.
% QCP   - Quadratic Complementarity Problem.
% MCP   - Mixed Complementarity Problem.
%
% formulated in AMPL.
%
% The interface detects the problem type, and sets the appropriate Prob
% structure.
% The AMPL problems have to be in the .nl format.
% The problem files could be generated in AMPL using the following commands:
%
% ampl -ogstub stub.mod stub.dat
% The .dat file is optional.
%
% amplAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% function Prob = amplAssign(nlFile, sparseFlag, solFile, HessPattern, x_L, x_U, x_0, ...
%                          setupFile, nProblem, pSepFunc, fLowBnd, VarWeight, KNAPSACK, ...
%                          fIP, xIP, b_L, b_U, ConsPattern, c_L, c_U, ...
%                          x_min, x_max, f_opt, x_opt, PriLev);
%
% INPUT (Call with at least one parameter, the nl File)
%
% nlFile      The AMPL nl file.
% sparseFlag  If 1 the problem is sparse, otherwise dense.
% solFile     If 1 a .sol file containing the solution will be created.
%             The file is readable by AMPL.
% HessPattern This option is only available for NLP and MINLP problems.
%             n x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the Hessian and ones indicate values that might
%             be non-zero. If empty indicates estimation of all elements
%             HessPattern is used when estimating the Hessian numerically.
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
%             The AMPL lower bounds are no longer valid if set.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
%             The AMPL upper bounds are no longer valid if set.
% x_0         Starting values, default nx1 zero vector.
%             The AMPL starting values are no longer valid if set.
%
% Note:       The number n of the unknown variables x is set from the AMPL
%             nl file.
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default values.
%
% setupFile   This option is only available for LP, MILP, QP and MIQP
%             problems.
%             The (unique) name as a TOMLAB Init File. If nonempty amplAssign
%             will create a executable m-file with this name and the given
%             problem defined as the first problem in this file.
%             See lp_prob.m, mip_prob.m, qp_prob.m, and qp_prob.m the TOMLAB
%             predefined Init Files.
%             If empty, no Init File is created.
%
% nProblem    This option is only available for LP, MILP, QP and MIQP
%             problems.
%             Number of problems, or problem number, to define in the setupFile
%             Not used if setupFile is empty.
%
%             nProblem = 1 ==> File is created to make it easy to edit new
%             problems into the file. Text are included on how to add new
%             problems. The given problem is set as number 1.
%             If isempty(nProblem) same as nProblem=1.
%
%             length(nProblem) > 1 ==> A file suitable for large test sets
%             are setup, where the problem definition is read from mat-files.
%             Statements for problems nProblem(1) to nProblem(2) are defined.
%             The given input is assumed to be nProblem(1), and the
%             corresponding mat-file is created.
%
%             If nProblem > 1. Additional problems are assumed, and the only
%             thing done is to create a mat-file with the problem.
%
%             If isempty(setupFile), nProblem is not used
%
% pSepFunc    This option is only available for NLP problems.
%             Number of subfunctions defined, if the function f is a partially
%             separable function.
%             The function f (and g and H) must check on Prob.PartSep.index
%             if Prob.PartSep.index == 0, compute the full function
%             if Prob.PartSep.index == i > 0, compute the i:th subfunction
%             This feature is only implemented in the solver sTrustr.
%
% fLowBnd     This option is only available for LP, MILP, QP and NLP
%             problems.
%             A lower bound on the function value at optimum.
%             A good estimate is not critical. Use [] if not known at all.
%
% VarWeight   This option is only available for MILP, MIQP and MINLP
%             problems.
%             Weight for each variable in the variable selection phase.
%             A lower value gives higher priority. Setting
%             Prob.MIP.VarWeight = c; for knapsack problems improve
%             convergence.
%
% KNAPSACK    This option is only available for MILP problems.
%             True if a knapsack problem is to be solved,
%             then a knapsack heuristic is used.
%
% fIP         An upper bound on the IP value wanted. Makes it possible to
%             cut branches and avoid node computations.
%
% xIP         The x-values giving the fIP value.
%
% L I N E A R   C O N S T R A I N T S
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
%             The AMPL values are no longer valid if set.
%
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%             The AMPL values are no longer valid if set.
%
% N O N L I N E A R   C O N S T R A I N T S
% c_L         Lower bound vector in nonlinear constraints b_L <= c(x) <= b_U.
%             The AMPL values are no longer valid if set.
% c_U         Upper bound vector in nonlinear constraints b_L <= c(x) <= b_U.
%             The AMPL values are no longer valid if set.
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%
% A D D I T I O N A L   P A R A M E T E R S
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
% PriLev  Print level of amplAssign.
%
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving.

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2003-2007 by Tomlab Optimization Inc., USA. $Release: 6.0.0$
% Written Jul 8, 2003.    Last modified Sep 24, 2007.

function Prob = amplAssign(nlFile, sparseFlag, solFile, HessPattern, x_L, x_U, x_0, ...
    setupFile, nProblem, pSepFunc, fLowBnd, VarWeight, KNAPSACK, fIP, xIP, b_L, b_U, ConsPattern, c_L, c_U, ...
    x_min, x_max, f_opt, x_opt, PriLev)
if nargin < 25, PriLev=[];
  if nargin < 24, x_opt=[];
    if nargin < 23, f_opt=[];
      if nargin < 22, x_max=[];
        if nargin < 21, x_min=[];
          if nargin < 20, c_U=[];
            if nargin < 19, c_L=[];
              if nargin < 18, ConsPattern=[];
                if nargin < 17, b_U=[];
                  if nargin < 16, b_L=[];
                    if nargin < 15, xIP=[];
                      if nargin < 14, fIP=[];
                        if nargin < 13, KNAPSACK=[];
                          if nargin < 12, VarWeight=[];
                            if nargin < 11, fLowBnd=[];
                              if nargin < 10, pSepFunc=[];
                                if nargin < 9, nProblem=[];
                                  if nargin < 8, setupFile=[];
                                    if nargin < 7, x_0=[];
                                      if nargin < 6, x_U=[];
                                        if nargin < 5, x_L=[];
                                          if nargin < 4, HessPattern=[];
                                            if nargin < 3, solFile=[];
                                              if nargin < 2, sparseFlag=[];
                                                if nargin < 1, nlFile=[];

end, end, end, end, end,
end, end, end, end, end,
end, end, end, end, end,
end, end, end, end, end,
end, end, end, end, end

if(isempty(PriLev))
    PriLev = 0;
end

if(isempty(nlFile))
    error('Problem file not specified.');
end

% Sets the .nl extention for the file unless it has already been set.

if isempty(findstr(nlFile, '.nl'))
    nlFile = strcat(nlFile, '.nl');
end

if(exist(nlFile,'file') ~= 2)
    error('Problem file %s does not exist.', nlFile);
end

% Clear the mex files, so that a new problem can be initialized.
clear spamfunc amplfunc amplqp;

% Sets [Primal initial guesses, lower bounds on the primal vbls, upper
% bounds on the primal vbls, dual initial guess, lower bounds on
% constrained bodies, upper bounds on constrained bodies.

% Also checking license

AMPL.sparseFlag = sparseFlag;
if AMPL.sparseFlag == 1
    spamfunc('lic');
    [AMPL.x_0, AMPL.x_L, AMPL.x_U, AMPL.v, AMPL.c_L, AMPL.c_U, AMPL.objType, AMPL.objConst, AMPL.nbv, AMPL.niv, AMPL.n_con, AMPL.nlc, AMPL.lnc, AMPL.n_obj, AMPL.nlo, AMPL.nzo, AMPL.cv] = spamfunc(nlFile);
else
    amplfunc('lic');
    [AMPL.x_0, AMPL.x_L, AMPL.x_U, AMPL.v, AMPL.c_L, AMPL.c_U, AMPL.objType, AMPL.objConst, AMPL.nbv, AMPL.niv, AMPL.n_con, AMPL.nlc, AMPL.lnc, AMPL.n_obj, AMPL.nlo, AMPL.nzo, AMPL.cv] = amplfunc(nlFile);
end

% Number of variables.
AMPL.n = length(AMPL.x_0);


% Conditioning AMPL.x_0, so that they are in the allowed range.
AMPL.x_0 = max( AMPL.x_L, min( AMPL.x_0, AMPL.x_U) );


% Setting the solFile flag to indicate that solution file is needed.
AMPL.solFile = solFile;


% Setting the linear constraint matrix A.

try
    if AMPL.sparseFlag == 1
        Prob.A = spamfunc('init', 0, AMPL.x_0+eps, zeros(AMPL.n_con,1));
    else
        Prob.A = amplfunc('init', 0, AMPL.x_0+eps, zeros(AMPL.n_con,1));
    end
catch
    error('Starting point need to be modified to initialize problem');
end

% , AMPL.APart, AMPL.Arows, AMPL.Acols
% objType:   -1 for maximize, 1 for minimize.
% objConst:  Constant term in objective function.
% nbv:       Number of binary variables.
% niv:       Number of integer variables.
% n_con:     Total number of constraints.
% nlc:       Total number of nonlinear constraints.
% lnc:       Number of linear network constraints.
% n_obj:     Number of objectives.
% nlo:       Number of nonlinear objectives.
% nzo:       Number of nonzeros in objective gradients.
% Jac:       The Jacobian, for linear break out of linear constraints.

% Number of linear constraints.
AMPL.lc = AMPL.n_con - AMPL.nlc;

AMPL.lin =  [];	    % i: linear constraints
AMPL.nlin = [];     % i: nonlinear constraints

% Iterate through all equality and inequality constraints.
for i = 1:AMPL.nlc
    AMPL.nlin = [AMPL.nlin i];
end
for i = AMPL.nlc + 1:1:AMPL.nlc + AMPL.lc
    AMPL.lin = [AMPL.lin i];
end

% Separating the linear and nonlinear lower bounds.
bl       = AMPL.c_L(AMPL.lin,:);
AMPL.c_L = AMPL.c_L(AMPL.nlin,:);
AMPL.b_L = bl;

% Separating the linear and nonlinear upper bounds.
bu       = AMPL.c_U(AMPL.lin,:);
AMPL.c_U = AMPL.c_U(AMPL.nlin,:);
AMPL.b_U = bu;

% Setting the integer indices for linear, quadratic and mixed nlp problems.

if(AMPL.nlc == 0 & AMPL.nlo == 0 & (AMPL.nbv > 0 | AMPL.niv > 0 ))
    % This is for MILP problems only.
    % Setting the integer variables.
    AMPL.nivIndices = zeros(AMPL.n,1);

    for i = AMPL.n - AMPL.niv + 1:1:AMPL.n
        AMPL.nivIndices(i) = 1;
    end

    % Setting the binary variables and checking AMPL.x_0, must be <= 1 and >= 0
    for i = AMPL.n - AMPL.niv - AMPL.nbv + 1:1:AMPL.n - AMPL.niv
        AMPL.nivIndices(i) = 1;
        if AMPL.x_L(i) < 0
            AMPL.x_L(i) = 0;
        end
        if AMPL.x_U(i) > 1
            AMPL.x_U(i) = 1;
        end
    end
elseif((AMPL.nlc > 0 | AMPL.nlo > 0) & (AMPL.nbv > 0 | AMPL.niv > 0))
    % This is for MINLP and MIQP problems.
    % Setting the integer variables.
    AMPL.nivIndices = [];

    % Setting the binary variables and checking AMPL.x_0, must be <= 1 and >= 0
    for i = AMPL.n - AMPL.niv - AMPL.nbv + 1:1:AMPL.n - AMPL.niv
        AMPL.nivIndices = [AMPL.nivIndices i];
        if AMPL.x_L(i) < 0
            AMPL.x_L(i) = 0;
        end
        if AMPL.x_U(i) > 1
            AMPL.x_U(i) = 1;
        end
    end

    for i = AMPL.n - AMPL.niv + 1:1:AMPL.n
        AMPL.nivIndices = [AMPL.nivIndices i];
    end
end

% Problem identification/classification.
% LP, MILP, QP, MIQP, NLP, MINLP

if (AMPL.nlc == 0 & AMPL.nlo == 0 & AMPL.nbv == 0 & AMPL.niv == 0)
    AMPL.problemType = 'LP';
elseif(AMPL.nlc == 0 & AMPL.nlo == 0 & (AMPL.nbv > 0 | AMPL.niv > 0 ))
    AMPL.problemType = 'MILP';
elseif(AMPL.nlc > 0 & AMPL.nbv == 0 & AMPL.niv == 0 )
    AMPL.problemType = 'NLP';
elseif(AMPL.nlc > 0 & (AMPL.nbv > 0 | AMPL.niv > 0 ))
    AMPL.problemType = 'MINLP';
elseif(AMPL.nlc == 0 & AMPL.nlo > 0 & AMPL.nbv == 0 & AMPL.niv == 0)
    % The problem is either QP or NLP.
    % Do QP check.
    amplqp('lic');
    qpcheck = amplqp(nlFile);
    clear amplqp;
    if qpcheck > 0
        AMPL.problemType = 'QP';
    else
        AMPL.problemType = 'NLP';
    end
elseif(AMPL.nlc == 0 & AMPL.nlo > 0 & (AMPL.nbv > 0 | AMPL.niv > 0))
    % The problem is either MIQP or MINLP
    % Do QP Check.
    amplqp('lic');
    qpcheck = amplqp(nlFile);
    clear amplqp;
    if qpcheck > 0
        AMPL.problemType = 'MIQP';
    else
        AMPL.problemType = 'MINLP';
    end
end

if any(AMPL.cv > 0)
    if any(strcmpi(AMPL.problemType,{'LP','MILP'}))
        AMPL.problemType = 'LCP';
    elseif any(strcmpi(AMPL.problemType,{'QP','MIQP'}))
        AMPL.problemType = 'QCP';
    else
        AMPL.problemType = 'MCP';
    end
end

Prob.AMPL = AMPL;

% Name of TOMLAB functions evaluating AMPL functions.
ampl_c   = 'ampl_c';
ampl_dc  = 'ampl_dc';
ampl_d2c = 'ampl_d2c';

% If no constraints, then no files should be named.
if isempty(AMPL.c_L) & isempty(AMPL.c_U)
    ampl_c   = [];
    ampl_dc  = [];
    ampl_d2c = [];
    ConsPattern = [];
end

switch AMPL.problemType

    % Assigning the problem to correct type.

    case 'LP'
        % Assigning LP.
        % Setting the objective function, c for a linear programming problem is
        % the derivative of f.
        Prob.c = ampl_g(AMPL.x_0,Prob);

        Prob = lpAssign(Prob.c, Prob.A, AMPL.b_L, AMPL.b_U, AMPL.x_L, AMPL.x_U, AMPL.x_0, nlFile,...
            setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt);

        if(PriLev > 0)
            fprintf('\n LP Assigned. \n');
        end

    case 'MILP'
        % Assigning MILP
        Prob.c = ampl_g(AMPL.x_0,Prob);

        Prob = mipAssign(Prob.c, Prob.A, AMPL.b_L, AMPL.b_U, AMPL.x_L, AMPL.x_U, AMPL.x_0, nlFile,...
            setupFile, nProblem, ...
            AMPL.nivIndices, VarWeight, KNAPSACK, fIP, xIP, ...
            fLowBnd, x_min, x_max, f_opt, x_opt);

        if(PriLev > 0)
            fprintf('\n MILP Assigned. \n');
        end

    case 'QP'
        % Assigning QP
        Prob.F = ampl_H(zeros(AMPL.n,1),Prob);
        Prob.c = ampl_g(zeros(AMPL.n,1),Prob);

        Prob = qpAssign(Prob.F, Prob.c, Prob.A, AMPL.b_L, AMPL.b_U, AMPL.x_L, AMPL.x_U, AMPL.x_0, nlFile,...
            setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt);

        if(PriLev > 0)
            fprintf('\n QP Assigned. \n');
        end

    case 'MIQP'
        %Assigning MIQP
        Prob.F = ampl_H(zeros(AMPL.n,1),Prob);
        Prob.c = ampl_g(zeros(AMPL.n,1),Prob);

        Prob = miqpAssign(Prob.F, Prob.c, Prob.A, AMPL.b_L, AMPL.b_U, AMPL.x_L, AMPL.x_U, AMPL.x_0,...
            AMPL.nivIndices, VarWeight, fIP, xIP,...
            nlFile, setupFile, nProblem, x_min, x_max, f_opt, x_opt);

        if(PriLev > 0)
            fprintf('\n MIQP Assigned. \n');
        end

    case 'NLP'
        % Assigning NLP
        Prob = conAssign('ampl_f', 'ampl_g', 'ampl_H', HessPattern, AMPL.x_L, AMPL.x_U, nlFile, AMPL.x_0, ...
            pSepFunc, fLowBnd, ...
            Prob.A, AMPL.b_L, AMPL.b_U, ampl_c, ampl_dc, ampl_d2c, ConsPattern, AMPL.c_L, AMPL.c_U, ...
            x_min, x_max, f_opt, x_opt);
        if ~isempty(c_L)
            Prob.c_L = c_L;
        end
        if ~isempty(c_U)
            Prob.c_U = c_U;
        end

        if(PriLev > 0)
            fprintf('\n NLP Assigned. \n');
        end

    case 'MINLP'
        % Assigning MINLP
        Prob = minlpAssign('ampl_f', 'ampl_g', 'ampl_H', HessPattern, AMPL.x_L, AMPL.x_U, nlFile, AMPL.x_0, ...
            AMPL.nivIndices, VarWeight, fIP, xIP,  ...
            Prob.A, AMPL.b_L, AMPL.b_U, ampl_c, ampl_dc, ampl_d2c, ConsPattern, AMPL.c_L, AMPL.c_U, ...
            x_min, x_max, f_opt, x_opt);
        if ~isempty(c_L)
            Prob.c_L = c_L;
        end
        if ~isempty(c_U)
            Prob.c_U = c_U;
        end

        if(PriLev > 0)
            fprintf('\n MINLP Assigned. \n');
        end

    case 'LCP'
        Prob.c = ampl_g(AMPL.x_0,Prob);
        numCC = nnz(AMPL.cv);
        MPEC = spalloc(numCC,6,numCC*2);
        count = 1;
        for i=1:AMPL.n_con
            if AMPL.cv(i) > 0
                MPEC(count,[1, 3]) = [AMPL.cv(i) i];
                count = count+1;
            end
        end

        Prob = lcpAssign(AMPL.objType*Prob.c, AMPL.x_L, AMPL.x_U, AMPL.x_0, Prob.A,AMPL.b_L, AMPL.b_U,...
            MPEC, nlFile, x_min, x_max);

        if(PriLev > 0)
            fprintf('\n LCP Assigned. \n');
        end
    case 'QCP'
        Prob.F = ampl_H(zeros(AMPL.n,1),Prob);
        Prob.c = ampl_g(zeros(AMPL.n,1),Prob);
        numCC = nnz(AMPL.cv);
        MPEC = spalloc(numCC,6,numCC*2);
        count = 1;
        for i=1:AMPL.n_con
            if AMPL.cv(i) > 0
                MPEC(count,[1, 3]) = [AMPL.cv(i) i];
                count = count+1;
            end
        end
        Prob = qcpAssign(Prob.F, Prob.c, Prob.A, AMPL.b_L, AMPL.b_U, AMPL.x_L, AMPL.x_U, AMPL.x_0,...
            MPEC, nlFile, x_min, x_max);

        if(PriLev > 0)
            fprintf('\n QCP Assigned. \n');
        end
    case 'MCP'
        numCC = nnz(AMPL.cv);
        MPEC = spalloc(numCC,6,numCC*2);
        count = 1;
        for i=1:AMPL.n_con
            if AMPL.cv(i) > 0
                if AMPL.cv(i) <= AMPL.lc
                    MPEC(count,[1, 3]) = [AMPL.cv(i) i];
                else
                    MPEC(count,[1, 5]) = [AMPL.cv(i) i-AMPL.lc];
                end
                count = count+1;
            end
        end
        Prob = mcpAssign('ampl_f', 'ampl_g', 'ampl_H', HessPattern, AMPL.x_L, AMPL.x_U, nlFile, AMPL.x_0, ...
               MPEC, fLowBnd, Prob.A, AMPL.b_L, AMPL.b_U, ampl_c, ampl_dc, ampl_d2c, ConsPattern, AMPL.c_L, AMPL.c_U, ...
               x_min, x_max, f_opt, x_opt);
        if ~isempty(c_L)
            Prob.c_L = c_L;
        end
        if ~isempty(c_U)
            Prob.c_U = c_U;
        end

        if(PriLev > 0)
            fprintf('\n MCP Assigned. \n');
        end

end % Switch

% Overriding AMPL values.
if ~isempty(x_L)
    Prob.x_L = x_L;
end
if ~isempty(x_U)
    Prob.x_U = x_U;
end
if ~isempty(x_0)
    Prob.x_0 = x_0;
end
if ~isempty(b_L)
    Prob.b_L = b_L;
end
if ~isempty(b_U)
    Prob.b_U = b_U;
end

% Setting the AMPL problem names, and Presolved variables.

% .col file, contains the variable names.
AMPL.vblNames = [];
colFile = strrep(nlFile, '.nl', '.col');
if(exist(colFile,'file') == 2)
    colFid = fopen(colFile);
    i = 1;
    tLine = fgetl(colFid);
    while tLine ~= -1
        AMPL.vblNames{i} = tLine;
        tLine = fgetl(colFid);
        i = i + 1;
    end
    fclose(colFid);
end

% .row file, contains the names of the constaints.
AMPL.nlconstrNames = [];
AMPL.lconstrNames = [];
rowFile = strrep(nlFile, '.nl', '.row');
if(exist(rowFile,'file') == 2)
    rowFid = fopen(rowFile);
    i = 1;
    tLine = fgetl(rowFid);
    while tLine ~= -1
        AMPL.nlconstrNames{i} = tLine;
        tLine = fgetl(rowFid);
        i = i + 1;
    end
    fclose(rowFid);
    % Setting the names of the objectives.
    AMPL.objNames = [];
    objNamesStart = length(AMPL.nlconstrNames) - AMPL.n_obj + 1;
    objNamesEnd = length(AMPL.nlconstrNames);
    AMPL.objNames = AMPL.nlconstrNames(objNamesStart:objNamesEnd);
    AMPL.nlconstrNames(objNamesStart:objNamesEnd) = [];
    AMPL.lconstrNames = AMPL.nlconstrNames(AMPL.nlc + 1:AMPL.nlc + AMPL.lc);
    AMPL.nlconstrNames(AMPL.lin) = [];
end

% .slc file, contains the names of the eliminated constraints.
AMPL.elimConstr = [];
slcFile = strrep(nlFile, '.nl', '.slc');
if(exist(slcFile,'file') == 2)
    slcFid = fopen(slcFile);
    i = 1;
    tLine = fgetl(slcFid);
    while tLine ~= -1
        AMPL.elimConstr{i} = tLine;
        tLine = fgetl(slcFid);
        i = i + 1;
    end
    fclose(slcFid);
end
AMPL = rmfield(AMPL, {'n', 'niv', 'nbv', 'nlc', 'lnc', 'n_obj', 'nlo', 'lc', 'lin', 'x_0', 'x_L', 'x_U', 'v', 'c_L', 'c_U', 'b_L', 'b_U'});
Prob.AMPL = AMPL;
% Reformulated problems that have passed through BuildMPEC need the AMPL
% structure in Prob.orgProb as well. 
if isfield(Prob,'orgProb')
   Prob.orgProb.AMPL = AMPL;
end

% MODIFICATION LOG
%
% 030708 medv Written.
% 030731 ango Fixed error messages and x_0 safeguard
% 030731 ango Changed local variables b_L,b_U to bl,bu to avoid
%             conflict with input arguments with same names
% 040519 ango Removed strange %d printouts
% 041111 med  No constraint files given if no constraints
% 050117 med  mlint review
% 050426 med  eps added to linear constraint initialization, error added
% 050426 med  objType added for LP, MILP, QP and MIQP
% 070901 med  LCP, QCP, MCP added
% 070924 ango Put AMPL structure in Prob.orgProb as well
