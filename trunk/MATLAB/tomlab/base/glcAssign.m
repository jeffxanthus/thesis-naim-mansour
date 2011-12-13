% glcAssign is a direct way of setting up a global mixed-integer programming
% (glc) problem in the TOMLAB Quick (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = glcAssign(...)
%
% It is then possible to solve the glc problem using a TOMLAB glc solver
% e.g.
%    Result = tomRun('glcSolve',Prob,2);
%    Result = tomRun('glcFast',Prob,2);
%    Result = tomRun('glcCluster',Prob,2);
%
% -----------------------------------------------------
%
% glc minimization problem:
%
%
%        min    f(x),  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%              c_L <= c(x) <= c_U
%
% Linear    equality equations: Set b_L==b_U
% Nonlinear equality equations: Set c_L==c_U
%
% Fixed     variables:          Set x_L==x_U. Both x_L and x_U must be finite.
%
% x(IntVars) are integer values, IntVars is an index set, a subset of [1:n].
%
% -----------------------------------------------------
%
% Syntax of glcAssign:
%
% function Prob = glcAssign(f, x_L, x_U, Name, A, b_L, b_U, ...
%                 c, c_L, c_U, x_0, IntVars, VarWeight, fIP, xIP, ...
%                 fLowBnd, x_min, x_max, f_opt, x_opt);
%
% f       Name of objective function f(x)
% x_L     Lower bounds on x, finite bounds must be given.
% x_U     Upper bounds on x, finite bounds must be given.
% Name    The name of the problem (string)
%
% The rest of the input parameters are optional
%
% A       The linear constraint matrix
% b_L     The lower bounds for the linear constraints, b_L <= A*x <= b_U.
% b_U     The upper bounds for the linear constraints, b_L <= A*x <= b_U.
% c       Name of constraint function c(x), computing nonlinear constraints
%         c(x) must be defined in a function, see e.g. glc_c.m
% c_L     Lower bound vector in nonlinear constraints, c_L <= c(x) <= c_U.
% c_U     Upper bound vector in nonlinear constraints, c_L <= c(x) <= c_U.
% x_0     Starting point x (may be empty, not used by DIRECT solvers, i.e.
%         glcFast, glcSolve, glbFast, glbSolve, glcCluster)
%
% b_L, b_U, c_L, c_U must either be empty or of full length
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
% Copyright (c) 2000-2007 by Tomlab Optimization Inc., Sweden. $Release: 6.0.0$
% Written Sept 28, 2000.  Last modified Sep 5, 2007.

function Prob = glcAssign(f, x_L, x_U, Name, A, b_L, b_U, ... 
                c, c_L, c_U, x_0, IntVars, VarWeight, fIP, xIP, ...
                fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 4
   error('glcAssign requires at least four parameters, f, x_L, x_U and Name'); 
end

if nargin < 20
   x_opt=[];
   if nargin < 19
      f_opt=[];
      if nargin < 18
         x_max=[];
         if nargin < 17
            x_min=[];
            if nargin < 16
               fLowBnd=[];
               if nargin < 15
                  xIP=[];
                  if nargin < 14
                     fIP=[];
                     if nargin < 13
                        VarWeight=[];
                        if nargin < 12
                           IntVars=[];
end, end, end, end, end, end, end, end, end

if nargin < 11
   x_0=[];
   if nargin < 10
      c_U=[];
      if nargin < 9
         c_L=[];
         if nargin < 8
            c=[];
            if nargin < 7
               b_U=[];
               if nargin < 6
                  b_L=[];
                  if nargin < 5
                     A=[];
end, end, end, end, end, end, end 

n  = max(length(x_L),length(x_U));

Prob          = ProbDef(1);
Prob.Name     = Name;
Prob.N        = n;
Prob.P        = 1;
Prob.probFile = 0;
Prob.PriLevOpt= 0;             % Set Print Level to 0 as default
Prob.x_opt    = x_opt;
Prob.f_opt    = f_opt;

if isempty(x_L) 
   error('x_L is empty! Only box-bounded problems allowed');
end

if isempty(x_U)
   error('x_U is empty! Only box-bounded problems allowed');
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

mA = size(A,1);
Prob.mLin = mA;

mN = max(length(c_L),length(c_U));
Prob.mNonLin = mN;

if mA+mN > 0
   Prob.probType = checkType('glc');
else
   Prob.probType = checkType('glb');
end

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
   fprintf('WARNING in glcAssign!!! ');
   fprintf('Constraint c is given. But no lower or upper bounds.\n');
end

if isempty(c) & (~isempty(c_L) | ~isempty(c_U))
   error('Constraint c not given, but lower and/or upper bounds.');
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

Prob = tomFiles(Prob, f, [],[], c);

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
% 000928 hkh  Created
% 011213 hkh  Changed Prob.MIP.SC into SC, SI, semi-continuous and semi-integer
% 011226 hkh  Changed Prob.MIP.SOS1 and SOS2 to sos1 and sos2
% 020822 hkh  Apostrophes missing in save statement
% 030210 hkh  Errors in comments about c_L and c_U
% 031201 hkh  Check if c defined and no bounds given, issue warning
% 040102 hkh  Add definition of fields mLin and mNonLin
% 040607 med  Help updates, extra checks removed
% 041123 hkh  Remove (weird) use of setupFile and nProblem
% 041123 hkh  Remove input KNAPSACK, set 0 in MIP structure
% 041201 hkh  Added check on columns of A, should be n
% 041201 hkh  There were two code parts doing roughly the same checks on A
% 041222 med  Checking lengths for x_L, x_U and x_0
% 050117 med  x_min and x_max were not properly set
% 050117 med  Local variable probType removed
% 051014 med  Checks if c given, when c_L or c_U are
% 051216 med  Function handles allowed as inputs
% 060206 med  Added checks for crossover bounds
% 060818 hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822 med  All vectors set to full and double
% 061213 med  Moved most input checks to checkAssign
% 070222 hkh  Revised IntVars input format and definition of Prob.MIP
% 070905 med  Function handle enabled
