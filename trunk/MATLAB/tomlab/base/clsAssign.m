% clsAssign implements the TOMLAB (TQ) format for mixed-integer unconstrained
% and constrained nonlinear least squares problems.
% It is also suitable to define vector valued function problems with
% a corresponding Jacobian matrix as derivative.
% For example minimax problems are solved with infSolve after using clsAssign
% to define the problem.
% L1 fitting problems are solved with L1Solve after using clsAssign
% to define the problem.
%
% clsAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% -----------------------------------------------------
%
% Mixed-integer constrained nonlinear (weighted) least squares problem:
%
%
%        min    f(x) = 0.5 * (r' * r),  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% Linear    equality equations: Set b_L==b_U
% Nonlinear equality equations: Set c_L==c_U
% Fixed     variables:          Set x_L==x_U. Both x_L and x_U must be finite.
% x(IntVars) are integer values, IntVars is an index set, a subset of [1:n].
%
% -----------------------------------------------------
%
% Syntax of clsAssign:
%
% function Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
%                           y, t, weightType, weightY, SepAlg, fLowBnd, ...
%                           A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ...
%                           x_min, x_max, f_opt, x_opt, ...
%                           IntVars, VarWeight, fIP, xIP);
%
% INPUT (Call with at least eight parameters)
%
% r           Name of function that computes the residual
% J           Name of function that computes the Jacobian m x n - matrix
% JacPattern  m x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the Jacobian and ones indicate values that might
%             be non-zero. If empty indicates estimation of all elements
%             JacPattern is used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, JacPattern==[]
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as  a nx1 Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
% Note 1:     The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specify at least one of these with correct length,
%             then the others are given default values
%
% y           m x 1 vector with observations y(t) to be fitted
%
% Note 2:     In the nonlinear least squares computations, the length of y
%             is used, and if weightType == 1, the y values are used to
%             weight the residual. However, the user should still subtract
%             the observations y from the model in the residual computations
% Note 3:     If y(t) is not used, then provide y = zeros(m,1), because Tomlab
%             needs the length m to setup memory correct for some solvers
%
% -----------------------------------------------------------------------------
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
% -----------------------------------------------------------------------------
%
% t           m x 1 vector with time values
% weightType  Type of weighting
%             0 = No weights (default if [])
%             1 = weight with absolute value of observation vector y (y(i) == 0 is tested for)
%             2 = weight with user given weight vector or weight matrix given
%                 as input weightY (see below).  Either length m x 1 or m x m 
%             3 = A user given function defined in Prob.LS.weightY will be
%                 called to compute the weights, defined as for weightType=2
%                 The function, say "DefWht" should be defined as
%                 function wLS = DefWht(x,r,Prob), function wLS = DefWht(x,r) or
%                 function wLS = DefWht(x)
% weightY     Either a user given function when weightType == 3, or
%             a vector or matrix of weights, m x 1 or m x m (dense or sparse)
%             The weights should be positive and multiplicative
%             weightY.*r  will be applied, if weightY is a vector
%             weightY*r   will be applied, if weightY is a matrix
% SepAlg      Flag if to use separable nonlinear least squares
%
% fLowBnd     A lower bound on the function value at optimum. Default -1E300
%             A good estimate is not critical. Use [] if not known at all.
%
% L I N E A R   C O N S T R A I N T S
% A           Matrix A in linear constraints b_L<=A*x<=b_U. Dense or sparse.
% b_L         Lower bound vector in linear constraints, b_L<=A*x<=b_U.
% b_U         Upper bound vector in linear constraints, b_L<=A*x<=b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
% dc          Name of function that computes the constraint Jacobian mN x n
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
% c_L         Lower bound vector in nonlinear constraints, c_L<=c(x)<=c_U.
% c_U         Upper bound vector in nonlinear constraints, c_L<=c(x)<=c_U.
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
%
% A D D I T I O N A L MIP P A R A M E T E R S
%
% If solving a mixed-integer nonlinear least squares problems
%
% IntVars     The set of integer variables. Can be given in one of two ways:
%
%             1) a vector of indices, e.g. [1 2 5].  If [], no integer variables
%
%             2) a 0-1 vector of length <= n=length(x) where nonzero elements
%                indicate integer variables
%
% VarWeight   Priorities for each variable in the variable selection phase
%             A lower value gives higher priority.
%
% fIP         An upper bound on the IP value wanted. Makes it possible to
%             cut branches and avoid node computations.
% xIP         The x-values giving the fIP value.
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Oct 12, 2000.    Last modified July 22, 2011.

function Prob = clsAssign(r, J, JacPattern, x_L, x_U, Name, x_0, ...
                          y, t, weightType, weightY, SepAlg, fLowBnd, ...
                          A, b_L, b_U, c, dc, ConsPattern, c_L, c_U, ... 
                          x_min, x_max, f_opt, x_opt, ...
                          IntVars, VarWeight, fIP, xIP)

if nargin < 29
   xIP=[];
   if nargin < 28
      fIP=[];
      if nargin < 27 
         VarWeight=[];
         if nargin < 26
            IntVars=[];
if nargin < 25
   x_opt=[];
   if nargin < 24
      f_opt=[];
      if nargin < 23
         x_max=[];
         if nargin < 22
            x_min=[];
            if nargin < 21
               c_U=[];
               if nargin < 20
                  c_L=[];
                  if nargin < 19
                     ConsPattern=[];
                     if nargin < 18
                        dc=[];
                        if nargin < 17
                           c=[];
                           if nargin < 16
                              b_U=[];
                              if nargin < 15
                                 b_L=[];
                                 if nargin < 14
                                    A=[];
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

if nargin < 13
   fLowBnd=[];
   if nargin < 12
      SepAlg=[];
      if nargin < 11
         weightY=[];
         if nargin < 10
            weightType=[];
            if nargin < 9
               t=[];

end, end, end, end, end, 

n = max([length(x_L),length(x_U),length(x_0)]);

Prob       = ProbDef(1);
Prob.P     = 1;
Prob.N     = n;
Prob.Name  = deblank(Name);

if isempty(x_min)
   Prob.x_min = -1*ones(n,1);
else
   Prob.x_min = x_min;
end
if isempty(x_max)
   Prob.x_max = 1*ones(n,1);
else
   Prob.x_max = x_max;
end

Prob.f_opt = f_opt;
Prob.x_opt = x_opt;

Prob.LS.y  = full(double(y(:)));
if isempty(y)
   disp('clsAssign: WARNING - empty y, e.g. solver NLSSOL will not work');
end
Prob.LS.t  = full(double(t(:)));

if isempty(weightType) 
   Prob.LS.weightType=0;
else
   Prob.LS.weightType=weightType;
end
if isempty(SepAlg) 
   Prob.LS.SepAlg=0;
else
   Prob.LS.SepAlg=SepAlg;
end

Prob.LS.weightY   = weightY;
Prob.JacPattern   = JacPattern;
Prob.ConsPattern  = ConsPattern;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);

mA = size(A,1);
mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if ~isempty(IntVars)
   Prob.probType = checkType('minlp');
elseif mA+mN > 0
   Prob.probType = checkType('cls');
else
   Prob.probType = checkType('ls');
end

Prob.probFile=0;

if mN > 0
   if isempty(c_L) 
      Prob.c_L=-Inf*ones(mN,1);
   else
      Prob.c_L=full(double(c_L(:)));
   end
   if isempty(c_U) 
      Prob.c_U=Inf*ones(mN,1);
   else
      Prob.c_U=full(double(c_U(:)));
   end
   if any(Prob.c_L>Prob.c_U)
      error('c_L and c_U have crossover values');
   end
end
if mN == 0 & ~isempty(c)
   fprintf('WARNING in clsAssign!!! ');
   fprintf('Constraint c is given. But no lower or upper bounds.\n');
end

if isempty(c) & (~isempty(c_L) | ~isempty(c_U))
   error('Constraint c not given, but lower and/or upper bounds.');
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

if any(IntVars == 0) | all(IntVars == 1)
   Prob.MIP = struct('IntVars',logical(abs(IntVars)),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
else
   Prob.MIP = struct('IntVars',full(double(abs(IntVars(:)))),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
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

Prob = tomFiles(Prob, 'ls_f', 'ls_g', 'ls_H', c, dc, [], r, J);

% MODIFICATION LOG
%
% 001012  hkh  Written, based on probAssign
% 011204  hkh  Set Prob.f_Low = 0
% 020409  hkh  Improve comments for vector valued functions
% 020411  hkh  Add comments about L1
% 030210  hkh  Errors in comments about c_L and c_U
% 031201  hkh  Check if c defined and no bounds given, issue warning
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040330  hkh  Demand input y, 8 parameters. Display warning if y empty
% 040414  hkh  Expand to handle MIP parameters, mixed-integer NLLS
% 040607  med  Help updates
% 040728  med  tomFiles used instead
% 041201  hkh  Added check number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  mlint review
% 051014  med  Checks if c given, when c_L or c_U are
% 051216  med  Function handles allowed as inputs
% 060206  med  Added checks for crossover bounds
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060818  hkh  Remove incorrect setting Prob.f_Low = 0; (cause not just NLLS)
% 060822  med  All vectors set to full and double
% 061213  med  Moved most input checks to checkAssign
% 070222  hkh  Revised IntVars input format and definition of Prob.MIP
% 070905  med  Function handle enabled
% 091023  med  Help updated
% 110722  hkh  Incorrect help for weightType and weightY, expanded comments
% 110722  hkh  Prob.LS.y incorrectly set to weightY, overwriting the y vector
% 110722  hkh  Set Prob.LS.weightY = weightY
