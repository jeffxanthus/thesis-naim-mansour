% simAssign implements the TOMLAB (TQ) format for unconstrained
% and constrained nonlinear programming problems, where both
% the function and the constraints need to be computed at the
% same time. This is often the case for simulation problems.
% If the problem is mixed-integer, set the last four input parameters.
%
% simAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% function Prob = simAssign(fc, gdc, Hd2c, HessPattern, x_L, x_U, Name, x_0, ...
%                           fLowBnd, A, b_L, b_U, ConsPattern, c_L, c_U, ...
%                           x_min, x_max, f_opt, x_opt, ...
%                           IntVars, VarWeight, fIP, xIP);
%
% INPUT (Call with at least seven parameters)
%
% fc          Name of the function that computes the function value f(x)
%             and the mN-vector of constraints c(x)
% gdc         Name of the function that computes the n x 1 gradient vector
%             and the mN by n matrix dc with the constraint gradients
% Hd2c        Name of the function that computes the n x n Hessian matrix and
%             the second part of the Lagrangian function (normally not used)
% HessPattern n x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the Hessian and ones indicate values that might
%             be non-zero. If empty indicates estimation of all elements
%             HessPattern is used when estimating the Hessian numerically.
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
% Note:       The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specify at least one of these with correct length,
%             then the others are given default values
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% fLowBnd     A lower bound on the function value at optimum. Default -1E300
%             A good estimate is not critical. Use [] if not known at all.
%
% L I N E A R   C O N S T R A I N T S
% A           mA x n matrix A, linear constraints b_L<=A*x<=b_U. Dense or sparse
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c_L         Lower bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% c_U         Upper bound vector in nonlinear constraints c_L <= c(x) <= c_U.
% ConsPattern mN x n zero-one sparse or dense matrix, where 0 values indicate
%             zeros in the constraint Jacobian and ones indicate values that
%             might be non-zero. Used when estimating the Jacobian numerically
%             and for efficient linear algebra operations and memory storage
%             in the solvers.
%             Estimated before solve, if Prob.LargeScale==1, ConsPattern==[]
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

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2001-2007 by Tomlab Optimization Inc., Sweden. $Release: 6.0.0$
% Written May 24, 2003.    Last modified Sep 5, 2007.

function Prob = simAssign(fc, gdc, Hd2c, HessPattern, x_L, x_U, Name, x_0, ...
                          fLowBnd, A, b_L, b_U, ConsPattern, c_L, c_U, ... 
                          x_min, x_max, f_opt, x_opt, ...
                          IntVars, VarWeight, fIP, xIP)

if nargin < 23
   xIP=[];
   if nargin < 22
      fIP=[];
      if nargin < 21 
         VarWeight=[];
         if nargin < 20
            IntVars=[];
if nargin < 19
   x_opt=[];
   if nargin < 18
      f_opt=[];
      if nargin < 17
         x_max=[];
         if nargin < 16
            x_min=[];
            if nargin < 15
               c_U=[];
               if nargin < 14
                  c_L=[];
                  if nargin < 13
                     ConsPattern=[];
                     if nargin < 12
                        b_U=[];
                        if nargin < 11
                           b_L=[];
                           if nargin < 10
                              A=[];
                              if nargin < 9
                                 fLowBnd=[];
                                 if nargin < 8
                                    x_0=[];
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

n = max([length(x_L),length(x_U),length(x_0)]);

Prob       = ProbDef(1);
Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
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

Prob.HessPattern  = HessPattern;
Prob.ConsPattern  = ConsPattern;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

mA           = size(A,1);
mN           = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if ~isempty(IntVars)
   Prob.probType = checkType('minlp');
elseif mA+mN > 0
   Prob.probType = checkType('con');
else
   Prob.probType = checkType('uc');
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

if any(IntVars == 0) | all(IntVars == 1)
   Prob.MIP = struct('IntVars',logical(abs(IntVars)),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
else
   Prob.MIP = struct('IntVars',full(double(abs(IntVars(:)))),...
                  'VarWeight',VarWeight,'KNAPSACK',0, 'fIP',fIP, 'xIP',xIP,...
                  'PI',[], 'SC',[], 'SI',[], 'sos1',[], 'sos2',[]);
end

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

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

Prob.PartSep.pSepFunc= 0;

Prob = tomFiles(Prob, 'sim_f', 'sim_g', '', 'sim_c', 'sim_dc', '', [],[],[], fc, gdc);

Prob.simType = 1;
if isempty(gdc) 
   Prob.ConsDiff = 1; 
   Prob.NumDiff  = 1;
end

% MODIFICATION LOG
%
% 030524  hkh  Written, based on conAssign
% 030525  hkh  Add sim_g and sim_dc interfaces (calling sim_gdc)
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040414  hkh  Expand to handle MIP parameters, mixed-integer NLLS
% 040607  med  Help updates
% 040728  med  tomFiles used instead
% 041010  hkh  Semi colon missing on line Prob.probType = checkType('minlp');
% 041201  hkh  Add check on number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  mlint revision
% 051216  med  Function handles allowed as inputs
% 060206  med  Added checks for crossover bounds
% 060705  med  Updated help
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822  med  All vectors set to full and double
% 061213  med  Moved most input checks to checkAssign
% 070222  hkh  Revised IntVars input format and definition of Prob.MIP
% 070905  med  Function handle enabled