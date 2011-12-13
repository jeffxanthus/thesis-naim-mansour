% bmiAssign provides a direct way of setting up linear semidefinite
% programming (SDP) problems with linear or bilinear matrix inequality
% constraints (LMI or BMI) in the TOMLAB (TQ) format.
%
% The information is put into the Tomlab input problem structure Prob.
%
% Prob = bmiAssign( .... )
%
% It is then possible to solve the SDP problem using PENSDP or PENBMI
%
% Result = tomRun('pensdp',Prob);
%
%  OR
%
% Result = tomRun('penbmi',Prob);
%
% bmiAssign can also take a predefined PENSDP problem structure as input
% and create a problem structure to be used by TOMLAB.
%
% See the demonstration example in sdpDemo.m in tomlab\examples
%
% -----------------------------------------------------
%
% The LMI/BMI minimization problem is defined as
% 
%        min   0.5 * x' * F * x + c' * x.  c and x in R^n
%         x
%
%   subject to 
%
%      x_L <=   x  <= x_U     x_L,x_U in R^n
%      b_L <= A*x  <= b_U     b_L,b_U in R^ml, A in R^[mlxn]
% 
% LMI: Q0^i + sum(i=1:n) (x_j*Q^i_j) <= 0, i = 1, ..., m 
%
% BMI: Q0^i + sum(i=1:n) (x_j*Q^i_j) + ...
%            + sum(k=1:n) sum(l=1:n) (x(k)*x(l)*K^i_kl) <= 0, i=1,...,m
%
% For each matrix constraint i, Q0^i, Q^i_j (and K^i_kl if BMI) must
% have same quadratic size (or empty)
%
% Equality equations: Set b_L==b_U
% Fixed variables:    Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of bmiAssign:
%
% function Prob = bmiAssign(F, c, BMI, A, b_L, b_U, x_L, x_U, x_0,  
%                           Name, fLowBnd);
%
%
% INPUT (One parameter F must always be given. Empty gives default)
%
% F             The matrix F in 0.5 x' F x in the objective
%               function OR a structure in the PENSDP Structure
%               Format (PSF).
%
% c             The vector c in c'x in the objective function.
%
% BMI           Structure array describing the LMI/BMI constraints.
%               i = 1, ..., m matrix inequalities
%
%   BMI(i).Q{j} are cell-arrays of p x p matrices (or empty)
%   BMI(i).Qidx are column vectors with the indices of the
%               corresponding matrix in BMI(i).Q (or empty)
%
%   BMI(i).K{j} are cell-arrays of p x p matrices (or empty)
%   BMI(i).Kidx are matrices with the indices of the corresponding
%               matrix in BMI(i).K. These matrices has two columns,
%               one for each index (or empty)
%
%               Each cell BMI(i).Q and BMI(i).K must be of the same
%               size p x p, but two inequalitíes BMI(i) and BMI(j)
%               may have differentyle sized matrices. Constraint i
%               is:
%
%   sum(j=1..k_L) (x(BMI(i).Qidx(j)) * BMI(i).Q{j}) +
%   + sum(j=1..k_B) (x(BMI(i).Kidx(j, 1)) * x(BMI(i).Kidx(j, 2)) *
%   * BMI(i).K{j}) <= 0
%   
%               where k_L is the number of non-zero LMI matrices
%               and k_B is the number of non-zero BMI
%               matrices. x(0) is treated as a constant 1 and is
%               not a variable, thus the constant matrix of the
%               LMI/BMI constraint is indexed with a zero (0) in 
%               the BMI(i).Qidx column vector.
%               
%               Any empty Q or K are interpreted as being all-zero.
%
% A             The linear constraint matrix 
% b_L           The lower bounds for the linear constraints
% b_U           The upper bounds for the linear constraints
%
% x_L           Lower bounds on x
% x_U           Upper bounds on x
%
%               b_L, b_U, x_L, x_U must either be empty or of full length
%
% x_0           Starting point x (may be empty)
%
% Name          The name of the problem (string)
%
% fLowBnd       A lower bound on the function value at optimum. Default -1E300
%               Use [] if not known at all. BMI solvers do not use fLowBnd

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written July 2, 1999.  Last modified Dec 13, 2006.

function Prob = bmiAssign(F, c, BMI, A, b_L, b_U, x_L, x_U, x_0, ...
			  Name, fLowBnd)

if nargin < 11
 fLowBnd = [];
 if nargin < 10
   Name=[];
   if nargin < 9
     x_0=[];
     if nargin < 8
       x_U=[];
       if nargin < 7
         x_L=[];
         if nargin < 6
           b_U=[];
           if nargin < 5
             b_L=[];
             if nargin < 4
               A=[];
               if nargin < 3
                 BMI=[];
                 if nargin < 2
                   c=[];
                   if nargin < 1
                     error('bmiAssign requires at least one parameter F'); 
end,end,end,end,end,end,end,end,end,end,end

probType = checkType('bmi');

Prob           =ProbDef(1);
Prob.P        = 1;
Prob.probFile = 0;
Prob.probType = probType;
Prob.Name     = deblank(Name);

if isstruct(F)
   Prob.x_0        = F.x0(:);
   Prob.N          = F.vars;
   Prob.QP.c       = F.fobj;
   Prob.QP.F       = [];
   Prob.PENOPT.pen = F;
   return
end

Prob.QP.F = F;
Prob.QP.c = c(:);
Prob.P    = 1;
n         = max(size(F,1),length(c));
Prob.N    = n;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
Prob.mNonLin = 0;

% Store in Prob.PENOPT.LMI
%Prob.PENOPT.LMI = LMI;

% Bi-linear parts of inequalities. No check at this point.
Prob.PENOPT.SDP = BMI;

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

Prob = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');

% MODIFICATION LOG
%
% 030127 ango Written (based on sdpAssign)
% 030524 hkh  Use tomFiles to define Prob.FUNCS structure
% 040102 hkh  Add definition of fields mLin and mNonLin
% 040607 med  Help updates
% 040728 med  tomFiles used instead
% 041201 hkh  Added check on number of columns in A, should be n
% 041222 med  Checking lengths for x_L, x_U and x_0
% 041214 frhe Changed to handle the new LMI/BMI format
% 041222 frhe Fixed bug checking for wrong structure input parameter.
% 041229 frhe Fixed bug in the Prob.PENOPT.SDP assignment.
% 050111 frhe Removed old LMI check.
% 050117 med  mlint review
% 050127 frhe nProblem and setupFile removed from argument list
% 060206 med  Added checks for crossover bounds
% 060818 hkh  Incorrect comments about fLowBnd not used, changed
% 060818 hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822 med  All vectors set to full and double
% 061213 med  Moved most input checks to checkAssign