% cgoAssign is a direct way of setting up a 
% costly (CPU-expensive) global mixed-integer black-box programming
% (cgo) problem in the TOMLAB Quick (TQ) format.
%
% If there are nonlinear constraints, the constraints might either be:
% 1. Costly constraints - computed at the same time as f(x). Often the
%    case for simulation problems.
% 2. Noncostly constraints - computed separately and treated by subsolvers
%    at the same time as considering linear constraints and box-bounds.
%    Derivatives for these noncostly constraints could be defined as
%    well as the pattern of the constraint Jacobian (ConsPattern).
% If the problem is mixed-integer or pure Integer-valued, IntVars, and 
% optionally  VarWeigtht, fIP and xIP are defined.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = cgoAssign(...)
%
% It is then possible to solve the cgo problem using a TOMLAB /CCO solver
% e.g.
%    Result = tomRun('rbfSolve',Prob,2);
%    Result = tomRun('ego',Prob,2);
%    Result = tomRun('arbfmip',Prob,2);
%
%
% -----------------------------------------------------
%
% CGO minimization problem:
%
%
%        min    f(x),  x in R^n
%         x
%        s/t   x_L  <=   x   <= x_U
%              b_L  <= A x   <= b_U
%              c_L  <= c(x)  <= c_U
%              Cc_L <= Cc(x) <= Cc_U
%
% Linear    equality equations: Set b_L(.) == b_U(.) for these.
% Nonlinear equality equations: Set c_L(.) == c_U(.) for these.
% Costly nonlinear equality constraints: Set Cc_L(.) == Cc_U(.) for these.
%
% Both x_L and x_U must be finite.
% Fixed     variables:   Set x_L(.) == x_U(.). 
%
% x(IntVars) are integer valued, IntVars is an index set, a subset of [1:n].
%
% -----------------------------------------------------
%
% Syntax of cgoAssign:
%
% function Prob = cgoAssign(fc, x_L, x_U, Name, Cc_L, Cc_U, A, b_L, b_U, ...
%                           c, dc, d2c, ConsPattern, c_L, c_U, x_0, ...
%                           IntVars, VarWeight, fIP, xIP, ...
%                           fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least four parameters)
%
% fc      Name of the function that computes the costly function value f(x)
%         and the costly mC-vector of constraints Cc(x)
% x_L     Lower bounds on x, finite bounds must be given.
% x_U     Upper bounds on x, finite bounds must be given.
% Name    The name of the problem (string)
%
% The rest of the input parameters are optional:
%
% Cc_L    Lower bound vector for costly nonlinear constraints Cc(x), 
% Cc_U    Upper bound vector for costly nonlinear constraints Cc(x), 
%         Cc_L <= Cc(x) <= Cc_U. 
%         If both Cc_L=[] & Cc_U=[], no costly constraints assumed, mC = 0;
%
% L I N E A R   C O N S T R A I N T S
% A       The linear constraint matrix
% b_L     The lower bounds for the linear constraints, b_L <= A*x <= b_U.
% b_U     The upper bounds for the linear constraints, b_L <= A*x <= b_U.
%
% N O N C O S T L Y   N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN noncostly nonlinear 
%             constraints c(x).
%             c(x) must be defined in a function, see e.g. glc_c.m
% dc          Name of function that computes the constraint Jacobian mN x n
% d2c         Name of function that computes the second part of the
%             Lagrangian function (only needed for some solvers)
%             See the help gateway routine nlp_d2c for an explanation of d2c
%
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
%
% c_L         Lower bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% c_U         Upper bound vector in nonlinear constraints c_L <= c(x) <= c_U.
%
% x_0         Starting point x (may be empty, not used by /CGO solvers)
%
% Cc_L, Cc_U, b_L, b_U, c_L, c_U must either be empty or of full length
%
% --------------------------------------------------------------------------
% The following variables are special for MIP problems, and are assigned to
% the field Prob.MIP.IntVars, etc.
%
% IntVars     The set of integer variables. Can be given in one of two ways:
%
%             1) a vector of indices, e.g. [1 2 5].  If [], no integer variables
%
%             2) a 0-1 vector of length <= n=length(x) where nonzero elements
%                indicate integer variables
%
% VarWeight  Weight for each variable in the variable selection phase.
%            A lower value gives higher priority. 
% fIP        An upper bound on the IP value wanted. Makes it possible to
%            cut branches and avoid node computations.
% xIP        The x-values giving the fIP value.
% --------------------------------------------------------------------------
% fLowBnd    A lower bound on the function value at optimum. Default -1E300
%            A good estimate is not critical. Use [] if not known at all.
% x_min      Lower bounds on each x-variable, used for plotting
% x_max      Upper bounds on each x-variable, used for plotting
% f_opt      Optimal function value(s), if known (Stationary points)
% x_opt      The x-values corresponding to the given f_opt, if known.
%            If only one f_opt, give x_opt as a 1 by n vector
%            If several f_opt values, give x_opt as a length(f_opt) by n matrix
%            If adding one extra column n+1 in x_opt,
%            0 indicates min, 1 saddle, 2 indicates max.
%            x_opt and f_opt is used in printouts and plots.

% Kenneth Holmstrom, Tomlab Optimization AB, E-mail: tomlab@tomopt.com
% Copyright (c) 2008-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written June 28, 2008.  Last modified June 28, 2008.

function Prob = cgoAssign(fc, x_L, x_U, Name, Cc_L, Cc_U, A, b_L, b_U, ...
                          c, dc, d2c, ConsPattern, c_L, c_U, x_0, ...
                          IntVars, VarWeight, fIP, xIP, ...
                          fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 4
   error('cgoAssign requires at least four parameters, f, x_L, x_U and Name'); 
end

if nargin < 25
   x_opt=[];
   if nargin < 24
      f_opt=[];
      if nargin < 23
         x_max=[];
         if nargin < 22
            x_min=[];
            if nargin < 21
               fLowBnd=[];
               if nargin < 20
                  xIP=[];
                  if nargin < 19
                     fIP=[];
                     if nargin < 18
                        VarWeight=[];
                        if nargin < 17
                           IntVars=[];
end, end, end, end, end, end, end, end, end

if nargin < 16
   x_0=[];
   if nargin < 15
      c_U=[];
      if nargin < 14
         c_L=[];
         if nargin < 13
            ConsPattern=[];
            if nargin < 12
               d2c=[];
               if nargin < 11
                  dc=[];
                  if nargin < 10
                     c=[];
                     if nargin < 9
                     b_U=[];
                        if nargin < 8
                        b_L=[];
                           if nargin < 7
                              A=[];
                              if nargin < 6
                                 Cc_U=[];
                                 if nargin < 5
                                    Cc_L=[];
end, end, end, end, end, end, end, end, end, end, end, end

if isempty(x_L) 
   error('x_L is empty! Only box-bounded problems allowed');
end

if isempty(x_U)
   error('x_U is empty! Only box-bounded problems allowed');
end

n             = max(length(x_L),length(x_U));
Prob          = ProbDef(1);
Prob.Name     = Name;
Prob.N        = n;
Prob.P        = 1;
Prob.probFile = 0;
Prob.PriLevOpt= 0;             % Set Print Level to 0 as default
Prob.x_opt    = x_opt;
Prob.f_opt    = f_opt;
Prob.ConsPattern  = ConsPattern;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);

if any(IntVars == 0) | all(IntVars == 1)
   Prob.MIP = struct('IntVars',logical(abs(IntVars)),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
else
   Prob.MIP = struct('IntVars',full(double(abs(IntVars(:)))),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
end

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end
Prob.PartSep.pSepFunc= 0;

mA           = size(A,1);
Prob.mLin    = mA;
mC           = max(length(Cc_L),length(Cc_U));
mN           = max(length(c_L),length(c_U));

if mA+mN+mC > 0
   Prob.probType = checkType('glc');
else
   Prob.probType = checkType('glb');
end

if mC > 0
   if isempty(Cc_L) 
      Cc_L=-Inf*ones(mC,1);
   else
      Cc_L=full(double(Cc_L(:)));
   end
   if isempty(Cc_U) 
      Cc_U=Inf*ones(mC,1);
   else
      Cc_U=full(double(Cc_U(:)));
   end
   if any(Cc_L > Cc_U)
      error('Cc_L and Cc_U have crossover values');
   end   
end
if mN > 0
   if isempty(c_L) 
      c_L=-Inf*ones(mN,1);
   else
      c_L=full(double(c_L(:)));
   end
   if isempty(c_U) 
      c_U=Inf*ones(mN,1);
   else
      c_U=full(double(c_U(:)));
   end
   if any(c_L > c_U)
      error('c_L and c_U have crossover values');
   end   
end
if mN == 0 & ~isempty(c)
   fprintf('WARNING in cgoAssign!!! ');
   fprintf('Constraint function c is given. But no lower or upper bounds.\n');
end

if isempty(c) & (~isempty(c_L) | ~isempty(c_U))
   error('Constraint function c not given, but lower and/or upper bounds.');
end

if mC == 0
   % Standard format
   Prob         = tomFiles(Prob, fc, [], [], c, dc, d2c);
   Prob.mNonLin = mN;
   if mN > 0
      Prob.c_L  = c_L;
      Prob.c_U  = c_U;
   end
elseif mN == 0
   % Pure simAssign format works out
   Prob = tomFiles(Prob,'sim_f','sim_g','','sim_c','sim_dc','', [],[],[],fc);
   Prob.simType = 1;
   Prob.mNonLin = mC;
   Prob.c_L     = Cc_L;
   Prob.c_U     = Cc_U;
else
   % mC > 0 & mN > 0
   % Both simAssign format needed, and also noncostly constraints
   Prob = tomFiles(Prob,'sim_f','sim_g','','sim_c','sim_dc','', [],[],[],fc);
   Prob.simType          = 2;
   Prob.c_L              = Cc_L;
   Prob.c_U              = Cc_U;
   % Store noncostly constraints
   Prob.CGO.c_L          = c_L;
   Prob.CGO.c_U          = c_U;
   Prob.CGO.c            = c;
   Prob.CGO.cNARG        = xnargin(c);
   Prob.CGO.dc           = dc;
   Prob.CGO.d2c          = d2c;
   Prob.CGO.ConsPattern  = ConsPattern;
   Prob.CGO.mNonLin      = mN;

   Prob.ConsPattern      = [];
   Prob.mNonLin          = mC;
end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print
if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

if isempty(x_min)
   Prob.x_min = -1*ones(n,1);
else
   Prob.x_min = x_min;
end
if isempty(x_max)
   Prob.x_max = ones(n,1);
else
   Prob.x_max = x_max;
end

% MODIFICATION LOG
%
% 080628 hkh  Created, from glcAssign, simAssign and conAssign
% 080710 hkh  Improve comments for VarWeight
