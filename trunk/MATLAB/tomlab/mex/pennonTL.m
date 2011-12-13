%  TOMLAB PENNON Solver
%
% function Result = pennonTL(Prob)
%
% -----------------------------------------------------------------------------
% INPUT:
% Prob   Problem structure in TOMLAB format
%
% Fields used in input structure Prob
% Use nlpsdpAssign to define the Prob structure,
% (or call Prob=ProbDef; and define each field)
%
%  x_L, x_U   Bounds on variables.
%  b_L, b_U   Bounds on linear constraints.
%  A          Linear constraint matrix.
%  c_L, c_U   Bounds on nonlinear constraints.
%
%  PriLevOpt  Print Level in pennonTL and MEX interface
%             (0 = off, 1 = summary, 2 = scalar information, 3 = verbose).
%
% -----------------------------------------------
% Fields used in Prob.PENOPT:
% -----------------------------------------------
%
%  NSDP         Structure array with Bilinear Matrix Inequality (BMI) data.
%
%    {i}.Y      p_i x p_i matrix variables (or empty)
%               Variable elements are indicated as NaN.
%               Constants with their respective values.
%    {i}.Y_L    p_i x p_i matrices (or empty)
%               representing the lower bounds on the elements in Y
%               If [] set as a p_i x p_i -Inf matrix.
%    {i}.Y_U    p_i x p_i matrices (or empty)
%               representing the upper bounds on the elements in Y
%               If [] set as a p_i x p_i Inf matrix.
%    {i}.lambda_L  The lower bound on the eigenvalue of NSDP(i).Y
%                  If [] set as a p_i x 1 zero vector.
%    {i}.lambda_U  The upper bound on the eigenvalue of NSDP(i).Y
%                  If [] set as a p_i x 1 Inf vector.
%    {i}.Type   Type of matrix constraint
%               0  Standard (Constraints may be infeasible)
%               1  Lower strict (The lower bounds will always be feasible)
%               2  Upper strict (The upper bounds will always be feasible)
%               3  Lower/upper strict (both bounds will always be
%                  feasible)
%               4  Slack (The block is a slack variable).
%
% Options   Options can be set in two ways.
%           Either setting the option arrays ioptions and doptions
%           in Prob.PENOPT.Options or by creating a field in Prob.PENOPT.Options:
%           Prob.PENOPT.Options.OPTIONNAME = OptionValue
%           The second method has precedence and is recommended.
%           Listed below is the placement of the options in the arrays
%           together with their respective field names.
%
% ioptions      (18x1)-vector with integer options.BIG
%               Set any element < 0 to make standard Tomlab parameter take
%               precedence.
%               Where applicable, standard Tomlab parameter is given in [].
%               Default values are PENNON defaults and are used
%               if nothing else is provided.
%
%  ioptions(1)  Maximum number of iterations of the overall algorithm.
%  MAXIT        [Prob.optParam.MaxIter] Default: 50
%
%  ioptions(2)  Maximum number of iterations of the inner loop.
%  NWITERS      [Prob.optParam.MinorIter] Default: 100
%
%  ioptions(3)  Print level: 0=silent, 1,2,3=summary,brief,full.
%  OUTLEV       Default: 0 (silent)
%
%  ioptions(4)  Hessian density check: 0=automatic check, 1=assume dense.
%  HESSIANMODE  Default: 0 (automatic)
%
%  ioptions(5)  Automatic scaling: 0=no, 1=yes
%  AUTOSCALE    Default: 0 (no)
%
%  ioptions(6)  Convex problem: 0=no, 1=yes
%  CONVEX       Default: 0 (no)
%
%  ioptions(7)  Treatment of equality constraints.
%  EQLTYMODE     0  As two inequalities, unsymmetric initialization
%                1  As two inequalities, symmetric initialization
%                2  Handled by standard augmented lagrangian
%                3  Direct handling (all equalities)
%               Default: 3 (direct)
%
%  ioptions(8)  Ignore inital solutions, 0=no, 1=yes
%  IGNOREINIT   Default: 0 (no)
%
%  ioptions(9)  Cholesky system mode
%  CHOLMODE      0  Solve directly
%                1  Solve augmented system
%                2  Split into two systems
%               Default: 0 (no)
%
%  ioptions(10) Newton stopping criteria
%  NWSTOPCRIT    0   || L(x^(k+1)) ||_2    < alfa
%                1   || L(x^(k+1)) ||_2    < alfa * || u_i^k - u_i^(k+1) ||_2
%                2   || L(x^(k+1)) ||_H^-1 < alfa * || L(x^k) ||_H^-1
%               Default: 2
%
%  ioptions(11) Penalty function
%  PENALTY       0  Logarithmic barrier and quadratic penalty
%                1  Reciprocal barrier and quadratic penalty
%               Default: 0 (logarithmic)
%
%  ioptions(12) Mode of solution of the Newton system.
%  NWTMODE       0  Cholesky method
%                1  Preconditioned conjugate gradient method
%                2  Preconditioned CG method with approximate
%                   Hessian calculation
%                3  Preconditioned CG method with user provided
%                   Hessian-vector routine
%                Default: 0 (Cholesky)
%
%  ioptions(13) Preconditioner type for the CG method.
%  PREC          0  No preconditioner
%                1  Diagonal
%                2  BFGS
%                3  Appoximate inverse
%                4  Symmetric Gauss-Seidel
%               Default: 0 (no preconditioner)
%
%  ioptions(14) Tuning parameter for Hessian assembling in Newton-mode 1-3
%  CMAXNZS      -1  Off
%               >0  On
%               Default: -1 (off)
%
%  ioptions(15) Automatic initialization of multipliers
%  AUTOINI       0  Off
%                1  Nonlinear (nonconvex) mode
%                2  LP/QP mode
%               Default: 0 (off)
%
%  ioptions(16) Perform penalty parameter update
%  USEPENUP      0  Adaptively
%                1  After each outer iteration
%               Default: 1 (each iteration)
%
%  ioptions(17) Box constraint mode
%  USEBARRIER    0  No special treatment
%                1  Use (strict) barrier function
%                2  Use (strict modified barrier function
%               Default: 0 (no special treatment)
%
%  ioptions(18) Derivative check
%  DERCHECK      0  No derivative check
%                1  Check gradients
%                2  Check Hessians
%               Default: 0 (no check)
%
% doptions      (14x1)-vector, floating point parameters.
%
%  doptions(1)  Stopping criterium for overall algorithm (rel. change in f).
%  PRECISION    [Prob.optParam.eps_f] Default: 1.0e-7
%
%  doptions(2)  Initial multiplier scaling factor.
%  UINIT        Default: 1.0
%
%  doptions(3)  Initial penalty value. Set lower (0.01-0.1) to maintain
%  PINIT        feasibility when starting from a feasible point.
%               Default: 1.0
%
%  doptions(4)  Stopping parameter alpha for the Newton/Trust region method
%  ALPHA        in the inner loop
%               Default 0.01
%
%  doptions(5)  Restriction factor of multiplier update
%  MU           Default: 0.5
%
%  doptions(6)  Penalty update
%  PENUP        Default: 0.1
%
%  doptions(7)  Minimal penalty
%  PEPS         Default: 1.0e-8
%
%  doptions(8)  Minimal multiplier
%  UMIN         Default: 1.0e-12
%
%  doptions(9)  Precision of the KKT conditions.
%  PRECKKT      Default: 1.0e-1
%
%  doptions(10) Minimum tolerance of the conjugate gradient algoritm
%  CGTOLMIN     Default: 5.0e-2
%
%  doptions(11) Update of tolerance of the conjugate gradient algorithm
%  CGTOLUP      Default: 1.0e0
%
%  doptions(12) Initial multiplier box constraints
%  UINITBOX     Default: 1.0e0
%
%  doptions(13) Initial multiplier nonlinear constraints
%  UINITNC      Default: 1.0e0
%
%  doptions(14) Initial multiplier sdp constraints
%  UINITSDP     Default: 1.0e0
%
% -----------------------------------------------------------------------------
% OUTPUT:
% Result  Structure with results (see ResultDef.m):
%
%   f_k       Function value at optimum or constraint deviation if infeasible.
%   x_k       Solution vector.
%   x_0       Initial solution vector.
%   g_k       Exact gradient computed at optimum.
%
%   xState    State of variables. Free == 0; On lower == 1; On upper == 2;
%             Fixed == 3;
%   bState    State of linear constraints. Free == 0; Lower == 1; Upper == 2;
%             Equality == 3;
%   v_k       Lagrangian multipliers (for bounds + dual solution vector).
%
%   ExitFlag  Exit status from pennon.m (similar to Tomlab).
%
%   Inform    PENNON information parameter.
%
%              0  - Solution obtained (No errors).
%              1  - Cholesky factorization of Hessian failed. The result may still be useful.
%              2  - No progress in objective value, problem probably infeasible.
%              3  - Linesearch failed. The result may still be useful.
%              4  - Maximum iteration limit exceeded. The result may still be useful.
%              5  - Wrong input parameters (ioptions,doptions).
%              6  - Memory error.
%              7  - Unknown error, please contact support@tomopt.com.
%
%   Iter      Number of iterations.
%   FuncEv    Number of function evaluations. Set to Iter.
%   GradEv    Number of gradient evaluations. Set to Iter.
%   ConstrEv  Number of constraint evaluations. Set to 0.
%   MinorIter Number of minor iterations. Always set to 0.
%
%   Solver    Name of the solver (PENNON).
%   SolverAlgorithm  Description of the solver.
%
% -----------------------------------------------
% Fields returned in Result.PENOPT:
%
%  NSDP         Structure array with solution matrices Y_k
%
%  iresults     Integer results:
%
%   iresults(1)  Number of outer iterations
%   iresults(2)  Number of inner iterations (Newton steps)
%   iresults(3)  Number of linesearch steps
%   iresults(4)  Elapsed time in seconds
%
%  fresults     Floating point results
%
%   fresults(1)  Objective f(x_k)
%   fresults(2)  Relative precision at x_k
%   fresults(3)  Feasibility at x_k
%   fresults(4)  Complementary slackness at x_k
%   fresults(5)  Gradient of augmented lagrangian at x_k
%
%  inform      PENNON inform value
%
% -----------------------------------------------------------------------
%
% For a problem description, see sdpnlpAssign.m
%
% -------------------------------------------------------------------------

% Bjorn Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2008 by Tomlab Optimization Inc., $Release: 5.0.0$
% Written Oct 2, 2003.  Last modified Oct 2, 2008.

function Result = pennonTL(Prob)

if nargin < 1, error('pennonTL needs the Prob structure as input');return;end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

% Need probType and solvType for PENNON
% Prob.solvType = 14; % BMI solver
% Prob = iniSolve(Prob,14,0,0);

Result=ResultDef(Prob);
Result.Solver='PENNON';
Result.SolverAlgorithm='NLP-SDP Solver PENNON 2.0';

x_0 = Prob.x_0;

if isfield(Prob.PENOPT,'pen')
   p = Prob.PENOPT.pen;
   if isempty(x_0)
      x_0 = p.x0(:);
   else
      p.x0 = x_0;
   end
else
   p = [];
end

% Define lower and upper bound arrays for PENNON
%
% Inf are changed to BIG (=1,0E38), -Inf to -BIG.
%
%   Used fields in structure Prob:
%     x_L      Lower bounds on x
%     x_U      Upper bounds on x
%     b_L      Lower bounds on linear constraints
%     b_U      Upper bounds on linear constraints
%     c_L      Lower bounds on nonlinear constraints
%     c_U      Upper bounds on nonlinear constraints

BIG=1.0E38;
[bl, bu, n, m, m2] = defblbu(Prob, BIG, 1);

if isempty(x_0)
   x_0 = zeros(n,1);
end

% Safe-guard for x_0
x_0  = max(bl(1:n,1),min(x_0(1:n,1),bu(1:n,1)));

% Set as constrained linear programming until any further classification.
% TODO: Check if derivative levels can be more precisely defined.
Prob = iniSolve(Prob,3,2,2);

% Check for struct Prob.PENOPT
Prob.PENOPT = DefPar(Prob, 'PENOPT');

% Defaults set in mex
ioptions = DefPar(Prob.PENOPT,'ioptions',[]);
doptions = DefPar(Prob.PENOPT,'doptions',[]);

[mA,nA] = size(Prob.A);

if ~isempty(Prob.A)
   % Relax this to allow linear constraints depending on matrix variables
   % and linear constraints defined in user callback functions
   %if nA~=n, error('Linear constraints A MUST have n columns!'); end
   %if mA~=m, error('Linear constraints A MUST have m rows!'); end
end

ixU = find(~isinf(Prob.x_U) & Prob.x_U <  BIG);
ixL = find(~isinf(Prob.x_L) & Prob.x_L > -BIG);

ibU = find(~isinf(Prob.b_U) & Prob.b_U <  BIG);
ibL = find(~isinf(Prob.b_L) & Prob.b_L > -BIG);

nU = length(ixU);
nL = length(ixL);

% total number of constraints
nconstr = m + m2;
nlin    = m;

% Bounds on variables. Standard (real) variables first,
% then vectorized matrix variables from Prob.PENOPT.NSDP.Y
lbv = Prob.x_L;
ubv = Prob.x_U;

% Bounds on constraints. Nonlinear before linear.
lbc = [Prob.c_L(:); Prob.b_L(:)];
ubc = [Prob.c_U(:); Prob.b_U(:)];

% No matrix variables unless defined in Prob.PENOPT.NSDP
nmv = 0;

% Matrix variables defined in cell array Prob.PENOPT.NSDP
NSDP = DefPar(Prob.PENOPT,'NSDP',[]);

% Initiate solver arguments determined from Prob.PENOPT.NSDP
nvars = n;
blks  = [];
lbmv  = [];
ubmv  = [];
mtype = [];
mnzs  = [];
mcol  = [];
mrow  = [];
xinit = x_0;
uinit = [];

n_Y = 0;

if isstruct(NSDP) 
   
   nmv = length(NSDP); % nmv instead of nsdp to avoid mixup with struct NSDP   
   
   for i=1:nmv
      
      % i as string for errors and warnings
      i_str = num2str(i);

      % Check for required field Y
      % Y defines the structure of matrix variable i
      if ~isfield(NSDP(i),'Y')
         error(strcat('Field Y missing in struct Prob.PENOPT.NSDP(',i_str,')'));
      else
         Y = NSDP(i).Y;
      end

      if isempty(Y)
         error(strcat('Field Y missing in struct Prob.PENOPT.NSDP(',i_str,')'));
      end
      % Y must be quadratic
      p = size(Y,1);
      if p ~= size(Y,2)
         error('Field Y in struct Prob.PENOPT.NSDP(',num2str(i),') must be quadratic');
      end

      % Y should be symmetric (but only upper triangle is extracted)
      Yt = Y';
      Y_num_ix = find(~isnan(Y));
      Y_NaN_ix = find(isnan(Y));  % Since NaN ~= NaN we need to compare indices
      Yt_NaN_ix = find(isnan(Yt));
      YT_NaN_ix = find(isnan(Y));
      if any(find(Y(Y_num_ix) ~= Yt(Y_num_ix)) | any(Y_NaN_ix ~= Yt_NaN_ix))
         warning(strcat('Field Y is not symmetric in struct Prob.PENOPT.NSDP(',i_str,')'));
      end

      % Find nonzero elements of matrix variables (both variables and nonzero constants)
      triu_Y = triu(Y);
      index_NaN = find(isnan(triu_Y));
      index_Nz = find(triu_Y ~= 0 & ~isnan(triu_Y));
      Y_var_index = sort([index_NaN; index_Nz]);
      
      n_Yi = length(Y_var_index);

      if n_Yi < 1
         error(strcat('Field Y in struct Prob.PENOPT.NSDP(',i_str,') is a zero matrix'));
      end    

      % Check Y_L and set -Inf if empty
      Y_L = DefPar(NSDP(i), 'Y_L');
      if isempty(Y_L)
         Y_L = Y;
         Y_L(index_NaN) = -Inf;
         Y_L(index_Nz) = Y(index_Nz);
      else
         if size(Y_L,1) ~= size(Y_L,2) | size(Y_L,1) ~= p
            error('Field Y_L in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y');
         end
      end
      
      % Check Y_U and set Inf if empty
      Y_U = DefPar(NSDP(i),'Y_U');
      if isempty(Y_U)
         Y_U = Y;
         Y_U(index_NaN) = Inf;
         Y_U(index_Nz) = Y(index_NaN);
      else
         if size(Y_U,1) ~= size(Y_U,2) | size(Y_U,1) ~= p
            error('Field Y_U in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y');
         end
      end
      
      % Check Y_0 and set 0 if empty
      Y_0 = DefPar(NSDP(i),'Y_0');
      if isempty(Y_0)
         Y_0 = Y;
         Y_0(Y_var_index) = min(max(0,Y_L(Y_var_index)),Y_U(Y_var_index));
      else
         if size(Y_0,1) ~= size(Y_0,2) | size(Y_0,1) ~= p
            error('Field Y_0 in struct Prob.PENOPT.NSDP(',i_str,') must be of the same size as Y');
         end
      end
      
      % Store variable bounds - Already done in nsdpAssign...
      %lbv = [lbv; Y_L(Y_var_index)];
      %ubv = [ubv; Y_U(Y_var_index)];

      % Store number of columns of each matrix variable
      blks(i) = p;

      % Lower and upper bounds on variable elements in matrix variables
      lbmv(i) = DefPar(NSDP(i),'lambda_L',-BIG);
      ubmv(i) = DefPar(NSDP(i),'lambda_U', BIG);

      % Type of matrix variables
      mtype(i) = DefPar(NSDP(i), 'Type', 0);

      % Number of nonzero elements in matrix variables
      mnzs(i) = length(Y_var_index);

      % Column and row indices of nonzero elements in matrix variables
      Y_index_matrix = zeros(blks(i));
      Y_index_matrix(Y_var_index) = 1;
      [ix, iy] = find(Y_index_matrix);
      mrow = [mrow; ix - 1]; % make indices zero-based
      mcol = [mcol; iy - 1]; % make indices zero-based

      % Update Prob.PENOPT.NSDP(i)
      Prob.PENOPT.NSDP(i).Y_L = Y_L;
      Prob.PENOPT.NSDP(i).Y_U = Y_U;
      Prob.PENOPT.NSDP(i).Y_0 = Y_0;
      Prob.PENOPT.NSDP(i).lambda_L = lbmv(i);
      Prob.PENOPT.NSDP(i).lambda_U = ubmv(i);
      Prob.PENOPT.NSDP(i).Type = mtype(i);

      % Indices will be reused in callbacks
      Prob.PENOPT.NSDP(i).Y_var_index = Y_var_index;
      
      % Count total number of variables in matrix variables
      n_Y = n_Y + n_Yi;
      
   end
else
   % No matrix variables
end

% Callback structures must be static
% Try to use ConsPattern and HessPattern
% If not defined, use full sets

%nnz_g = DefPar(Prob.PENOPT,'nnz_gradient',nvars);
nnz_g = nvars; % always dense for now

nnz_h = nnz(DefPar(Prob,'HessPattern',[]));
if nnz_h < 1
   [hp_ix, hp_iy] = find(tril(ones(nvars)));
   nnz_h = length(hp_ix);
else
   [hp_ix, hp_iy] = find(tril(Prob.HessPattern));
   nnz_h = length(hp_ix);
end
% zero-based indices for c-mex
hp_ix = hp_ix - 1;
hp_iy = hp_iy - 1;

nnz_dc = nnz(DefPar(Prob,'ConsPattern',[]));

% Constraints are returned one at a time
% Use a minimum static structure by eliminating any zero columns
% due to variables not involved in any constraints
if nnz_dc < 1
   nnz_dc = nvars*(nconstr-mA);
   [cp_ix, cp_iy] = find(ones(nconstr-mA,nvars));
   if nconstr-mA > 0
      cp_i = 1:nvars;
   else
      cp_i = [];
   end
else
   [cp_ix, cp_iy] = find(Prob.ConsPattern);
   [dummy, cp_i] = find(sum(Prob.ConsPattern,1));
end

% create a static d2cpattern based on the ConsPattern
cHp_ix = [];
cHp_iy = [];
cHp_sz = [1];

for i = 1:nconstr
   dcP_yi_i =  cp_iy(find(cp_ix == i));
   cHp_i = zeros(nvars);
   cHp_i(dcP_yi_i,dcP_yi_i) = 1;
   [cHp_ix_i,cHp_iy_i] = find(tril(cHp_i));
   cHp_ix = [cHp_ix; cHp_ix_i];
   cHp_iy = [cHp_iy; cHp_iy_i];
   cHp_sz(i+1) = cHp_sz(i) + length(cHp_ix_i);
end

% zero-based indices for c-mex
cp_i = cp_i - 1;
cp_ix = cp_ix - 1;
cp_iy = cp_iy - 1;

cHp_ix = cHp_ix - 1;
cHp_iy = cHp_iy - 1;
cHp_sz = cHp_sz - 1;

% add matrix variables to xinit - Already done in nsdpAssign
%xinit = [Prob.x_0; xinit];

% Callbacks need to separate standard variables from matrix variables
%Prob.PENOPT.n_X = n;
n_X = n - n_Y;

% uinit is not mentioned as a direct input vector in the manual,
% but exists as such in the Matlab-MEX from PENOPT.
% Any relation to the options uinit, uinitbox, uinitnc, uinitsdp?
% For now, leave empty, will result in a zero-vector in the MEX.

% Print level in MEX
PriLev=DefPar(Prob,'PriLevOpt',0);

% Adjust everything to be in the interval [-BIG, BIG]
ix=find(lbv < -BIG);  % x_L
lbv(ix)=-BIG;
ix=find(ubv > BIG);   % x_U
ubv(ix)=BIG;
ix=find(lbc < -BIG);  % b_L c_L
lbc(ix)=-BIG;
ix=find(ubc > BIG);   % b_U c_U
ubc(ix)=BIG;
ix=find(lbmv < -BIG); % lambda_L
lbmv(ix)=-BIG;
ix=find(ubmv > BIG);  % lambda_U
ubmv(ix)=BIG;


% Might be necessary to pass through standard callbacks (nlp_*.m):
% Update Prob with matrix variables merged into x.
Prob.N   = nvars;
Prob.x_0 = xinit;
Prob.x_L = lbv;
Prob.x_U = ubv;

uinit = [];

if m > 0
   % extend A to size m x nvars
   Prob.A = [Prob.A zeros(mA,max(nvars-nA,0))];
end

% Check for ioptions and doptions
Prob.PENOPT.Options = DefPar(Prob.PENOPT, 'Options', []);
ioptions = DefPar(Prob.PENOPT.Options, 'ioptions',-1*ones(18,1));
doptions = DefPar(Prob.PENOPT.Options, 'doptions',-1.0*ones(14,1));
ioptions = [ioptions(:); zeros(max(0,18-length(ioptions)),1) ];
doptions = [doptions(:); zeros(max(0,14-length(doptions)),1) ];

% Check for fields with options in Prob.PENOPT.Options
ioptions(1) =  DefPar(Prob.PENOPT.Options,'MAXIT',       ioptions(1));
ioptions(2) =  DefPar(Prob.PENOPT.Options,'NWITERS',     ioptions(2));
ioptions(3) =  DefPar(Prob.PENOPT.Options,'OUTLEV',      ioptions(3));
ioptions(4) =  DefPar(Prob.PENOPT.Options,'HESSIANMODE', ioptions(4));
ioptions(5) =  DefPar(Prob.PENOPT.Options,'AUTOSCALE',   ioptions(5));
ioptions(6) =  DefPar(Prob.PENOPT.Options,'CONVEX',      ioptions(6));
ioptions(7) =  DefPar(Prob.PENOPT.Options,'EQLTYMODE',   ioptions(7));
ioptions(8) =  DefPar(Prob.PENOPT.Options,'IGNOREINIT',  ioptions(8));
ioptions(9) =  DefPar(Prob.PENOPT.Options,'CHOLMODE',    ioptions(9));
ioptions(10) = DefPar(Prob.PENOPT.Options,'NWTSTOPCRIT', ioptions(10));
ioptions(11) = DefPar(Prob.PENOPT.Options,'PENALTY',     ioptions(11));
ioptions(12) = DefPar(Prob.PENOPT.Options,'NWTMODE',     ioptions(12));
ioptions(13) = DefPar(Prob.PENOPT.Options,'PREC',        ioptions(13));
ioptions(14) = DefPar(Prob.PENOPT.Options,'CMAXNZS',     ioptions(14));
ioptions(15) = DefPar(Prob.PENOPT.Options,'AUTOINI',     ioptions(15));
ioptions(16) = DefPar(Prob.PENOPT.Options,'PENUP',       ioptions(16));
ioptions(17) = DefPar(Prob.PENOPT.Options,'USEBARRIER',  ioptions(17));
ioptions(18) = DefPar(Prob.PENOPT.Options,'DERCHECK',    ioptions(18));

doptions(1) =  DefPar(Prob.PENOPT.Options,'PRECISION',   doptions(1));
doptions(2) =  DefPar(Prob.PENOPT.Options,'UINIT',       doptions(2));
doptions(3) =  DefPar(Prob.PENOPT.Options,'PINIT',       doptions(3));
doptions(4) =  DefPar(Prob.PENOPT.Options,'ALPHA',       doptions(4));
doptions(5) =  DefPar(Prob.PENOPT.Options,'MU',          doptions(5));
doptions(6) =  DefPar(Prob.PENOPT.Options,'PENUP',       doptions(6));
doptions(7) =  DefPar(Prob.PENOPT.Options,'PEPS',        doptions(7));
doptions(8) =  DefPar(Prob.PENOPT.Options,'UMIN',        doptions(8));
doptions(9) =  DefPar(Prob.PENOPT.Options,'PRECKKT',     doptions(9));
doptions(10) = DefPar(Prob.PENOPT.Options,'CGTOLMIN',    doptions(10));
doptions(11) = DefPar(Prob.PENOPT.Options,'CGTOLUP',     doptions(11));
doptions(12) = DefPar(Prob.PENOPT.Options,'UINITBOX',    doptions(12));
doptions(13) = DefPar(Prob.PENOPT.Options,'UINITNC',     doptions(13));
doptions(14) = DefPar(Prob.PENOPT.Options,'UINITSDP',    doptions(14));

% Some options are equivalent to standard TOMLAB parameters
if ioptions(1) < 0  ioptions(1) = Prob.optParam.MaxIter;   end
if ioptions(2) < 0  ioptions(2) = Prob.optParam.MinorIter; end

% THE CALL
[Obj, x_k, c_k, info, iresults, fresults] = ...
   Tpennon(nvars, nlin, nconstr, nmv, blks, nnz_g, nnz_h, lbv, ubv, ...
   lbc, ubc, lbmv, ubmv, mtype, mnzs, mrow, mcol, xinit, uinit, Prob.A, ...
   hp_ix, hp_iy, cp_i, cp_ix, cp_iy, cHp_sz, cHp_ix, cHp_iy,...
   ioptions, doptions, Prob.PriLevOpt, Prob);

Inform = info;
Iter = iresults(1);

switch Inform 
   case {2,3},        ExitFlag = 4;  % Infeasible
   case 4    ,        ExitFlag = 2;  % Unbounded
   case 5    ,        ExitFlag = 1;  % Too many iterations
   case {9,10,11,12}, ExitFlag = -1; % Fatal errors
   case 8    ,        ExitFlag = 10; % Input error
   otherwise ,        ExitFlag = 0;
end

switch Inform
   case 0 , Text = 'Solution obtained.';
   case 1 , Text = 'Suboptimal solution (large gradient)';
   case 2 , Text = 'Solution primal infeasible.';
   case 3 , Text = 'No progress, problem may be infeasible.';
   case 4 , Text = 'Primal unbounded or initial multipliers too small.';
   case 5 , Text = 'Maximum iteration limit exceeded. The result may still be useful.';
   case 6 , Text = 'Line search failed. The result may still be useful.';
   case 7 , Text = 'Cholesky factorization of Hessian failed. The result may still be useful.';
   case 8 , Text = 'Wrong input parameters.';
   case 9 , Text = 'Resource limit.';
   case 10, Text = 'Internal error, please contact support@tomopt.com';
   case 11 ,Text = 'Memory error, please contact support@tomopt.com.';
   case 12, Text = 'Unknown error, please contact support@tomopt.com.';
   otherwise Text = 'NOTE: UNKNOWN PENNON Inform value.';
end

Result.ExitText = Text;

Result.f_k = Obj;
Result.x_k = x_k(1:n);
Result.x_0 = x_0;
Result.c_k = c_k(1:m2);


Result.FuncEv    = Iter;
Result.GradEv    = Iter;
Result.ConstrEv  = Iter;
Result.Iter      = iresults(1);
Result.MinorIter = iresults(2);
Result.ExitFlag  = ExitFlag;
Result.Inform    = Inform;

Result.CPUtime   = iresults(4);
Result.REALtime  = iresults(4);

Result.PENOPT.iresults = iresults;
Result.PENOPT.fresults = fresults;
Result.PENOPT.info     = info;

x_pos = n_X;
for i = 1:nmv
   Y_i = Prob.PENOPT.NSDP(i).Y;
   Yi = Prob.PENOPT.NSDP(i).Y_var_index;
   nYi = length(Yi);

   % Extract upper triangle
   Y_i = triu(Y_i);

   % Update upper triangle
   Y_i(Yi) = x_k(x_pos+1:x_pos+nYi);

   % Make Y_k symmetric
   Result.PENOPT.NSDP(i).Y_k = Y_i + Y_i' - diag(diag(Y_i));

   x_pos = x_pos + nYi;
end

if mA > 0
   Ax = Prob.A * x_k;
else
   Ax = [];
end

% Compute Result.xState
Result = StateDef(Result, x_k, Ax, [],Prob.optParam.xTol, ...
   Prob.optParam.bTol, [], [lbv;lbc], [ubv;ubc],1);

if PriLev > 0
   fprintf('\n\n-->-->-->-->-->-->-->-->-->-->');
   fprintf('\nPENNON solving Problem %d:\n',Prob.P);
   fprintf('-->-->-->-->-->-->-->-->-->-->\n\n');
   fprintf('PENNON: Inform = %2d, ',Inform)
   fprintf('\n');
   for i = 1:size(Text,1)
      fprintf('%s',Text(i,:))
      fprintf('\n')
   end
   fprintf('\n');

   fprintf('Objective function at solution x %36.18f\n\n',Obj);
   fprintf('Iterations      %7d. ',Iter);
   fprintf('\n');

   if PriLev > 1
      if isempty(MAX_x)
         MAX_x=length(x_k);
      end
      fprintf('Optimal x = \n');
      xprinte(x_k(1:min(length(x_k),MAX_x)),'x_k:  ');
   end
end

Result=endSolve(Prob,Result);

% MODIFICATION LOG:
%
% 081002 bjo  Created
