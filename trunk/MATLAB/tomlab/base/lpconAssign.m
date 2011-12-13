% lpconAssign implements the TOMLAB (TQ) format for nonlinearly constrained
% linear programming problems.
%
% lpconAssign sets the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% This routine is usable if you want to call a solver directly and do not
% want to use the TOMLAB Init File format to define the problem.
%
% -----------------------------------------------------
%
% Nonlinearly constrained LP minimization problem:
%
%
%        min    d' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% Linear    equality equations: Set b_L==b_U
% Nonlinear equality equations: Set c_L==c_U
% Fixed     variables:          Set x_L==x_U. Both x_L and x_U must be finite.
%
% -----------------------------------------------------
%
% Syntax of lpconAssign:
%
% function Prob = lpconAssign(d, x_L, x_U, Name, x_0, A, b_L, b_U,...
%                             c, dc, d2c, ConsPattern, c_L, c_U,...
%                             fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least three parameters)
%
% d           The vector d in d'x in the objective function
%
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as a nx1  Inf vector.
%
% Note:       The number n of the unknown variables x are taken as
%             max(length(x_L),length(x_U),length(x_0))
%             You must specify at least one of these with correct length,
%             then the others are given default values
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
% L I N E A R   C O N S T R A I N T S
% A           mA x n matrix A, linear constraints b_L <= A*x <= b_U. Dense or sparse
% b_L         Lower bound vector in linear constraints b_L <= A*x <= b_U.
% b_U         Upper bound vector in linear constraints b_L <= A*x <= b_U.
%
% N O N L I N E A R   C O N S T R A I N T S
% c           Name of function that computes the mN nonlinear constraints
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
% A D D I T I O N A L   P A R A M E T E R S
% fLowBnd     A lower bound on the function value at optimum. Default -1E300
%             A good estimate is not critical. Use [] if not known at all.
%             Only used by some nonlinear solvers
%
% x_min   Lower bounds on each x-variable, used for plotting
% x_max   Upper bounds on each x-variable, used for plotting
% f_opt   Optimal function value(s), if known (Stationary points)
% x_opt   The x-values corresponding to the given f_opt, if known.
%         If only one f_opt, give x_opt as a 1 by n vector
%         If several f_opt values, give x_opt as a length(f_opt) by n matrix
%         If adding one extra column n+1 in x_opt, 0 is min, 1 saddle, 2 is max.
%         x_opt and f_opt is used in printouts and plots.
%
%
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2007 by Tomlab Optimization Inc., Sweden. $Release: 6.0.0$
% Written Nov 14, 2004.    Last modified Sep 5, 2007.

function Prob = lpconAssign(d, x_L, x_U, Name, x_0, ...
                A, b_L, b_U, c, dc, d2c, ConsPattern, c_L, c_U, ... 
                fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 19
   x_opt=[];
   if nargin < 18
      f_opt=[];
      if nargin < 17
         x_max=[];
         if nargin < 16
            x_min=[];
            if nargin < 15
               fLowBnd=[];
               if nargin < 14
                  c_U=[];
                  if nargin < 13
                     c_L=[];
                     if nargin < 12
                        ConsPattern=[];
                        if nargin < 11
                           d2c=[];
                           if nargin < 10
                              dc=[];
                              if nargin < 9
                                 c=[];
                                 if nargin < 8
                                    b_U=[];
                                    if nargin < 7
                                       b_L=[];
                                       if nargin < 6
                                          A=[];
                                          if nargin < 5
                                             x_0=[];
                                             if nargin < 4
                                                Name=[];
                                                if nargin < 3
                    error('lpconAssign requires at least three parameters');
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

n = max([length(x_L),length(x_U),length(x_0)]);
if ~isempty(A)
   if size(A,2) ~= length(d)
      error('Illegal length of d or number of columns in A.');
   end
end
Prob       = ProbDef(1);
Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);

Prob.QP.c = full(double(d(:)));
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

Prob.ConsPattern  = ConsPattern;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

mA = size(A,1);
mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if mA+mN > 0
   Prob.probType = checkType('con');
else
   Prob.probType = checkType('lp');
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
   fprintf('WARNING in conAssign!!! ');
   fprintf('Constraint c is given. But no lower or upper bounds.\n');
end

if isempty(c) & (~isempty(c_L) | ~isempty(c_U))
   error('Constraint c not given, but lower and/or upper bounds.');
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

Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H', c, dc, d2c);

% MODIFICATION LOG
%
% 041114  med  Written, based on conAssign
% 041201  hkh  Add check on number of columns of A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  mlint revision
% 050613  hkh  Change f to d for the coefficients in the objective
% 051014  med  Checks if c given, when c_L or c_U are
% 051212  med  More size checks added on inputs
% 051216  med  Function handles allowed as inputs
% 060206  med  Added checks for crossover bounds
% 060705  med  Updated help
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822  med  All vectors set to full and double
% 061213  med  Moved most input checks to checkAssign
% 070905  med  Function handle enabled