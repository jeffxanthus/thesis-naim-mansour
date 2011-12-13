% llsAssign implements the TOMLAB (TQ) format for the linear least
% squares problem.
%
% llsAssign is setting the variables normally needed for an optimization in
% the TOMLAB structure Prob.
%
% -----------------------------------------------------
%
% LLS minimization problem:
%
%
%        min   f(x) =  0.5 * ||Cx - y||^2.  x in R^n, || || L_2-norm
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of llsAssign:
%
% function Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
%                           t, weightType, weightY, ...
%                           A, b_L, b_U, ...
%                           x_min, x_max, f_opt, x_opt);
%
% INPUT (Call with at least five parameters)
%
% C           Matrix m x n in objective ||Cx - y||^2
% y           Vector m x 1 with observations in objective ||Cx - y||
% x_L         Lower bounds on parameters x. If [] set as a nx1 -Inf vector.
% x_U         Upper bounds on parameters x. If [] set as  a nx1 Inf vector.
% Name        The name of the problem (string)
% x_0         Starting values, default nx1 zero vector
%
%             The following parameters are optional, and problem type dependent
%             Set empty to get default value
%
% t           m x 1 vector with time values (if empty assumed to be 1:m)
%             C = C(t) and y = y(t)
% weightType  Type of weighting
%             0 = No weights (default if [])
%             1 = weight with absolute value of observation vector y (y(i) == 0 is tested for)
%             2 = weight with user given weight vector or weight matrix given
%                 as input weightY (see below)
%                 Must either have length m x 1 or m x m (could be sparse)
% weightY     Vector or matrix of weights, m x 1 or m x m (dense or sparse)
%             The weights should be positive and multiplicative
%             weightY*.y and diag(weightY)*C will be applied, if weightY is a vector
%             weightY*y and (weightY*C will be applied, if weightY is a matrix
% NOTE!       Any weighting is directly applied to C and y and the weighted results
%             are stored in Prob.LS.C and Prob.LS.y
%
% L I N E A R   C O N S T R A I N T S
% A           Matrix A in linear constraints b_L<=A*x<=b_U. Dense or sparse.
% b_L         Lower bound vector in linear constraints b_L<=A*x<=b_U.
% b_U         Upper bound vector in linear constraints b_L<=A*x<=b_U.
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
% Set the variable as empty if this variable is not needed for the particular
% kind of problem you are solving

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2002-2011 by Tomlab Optimization Inc., $Release: 7.8.0$
% Written Nov 8, 2000.    Last modified July 22, 2011.

function Prob = llsAssign(C, y, x_L, x_U, Name, x_0, ...
                          t, weightType, weightY, ...
                          A, b_L, b_U, x_min, x_max, f_opt, x_opt)

if nargin < 16
   x_opt=[];
   if nargin < 15
      f_opt=[];
      if nargin < 14
         x_max=[];
         if nargin < 13
            x_min=[];
            if nargin < 12
               b_U=[];
               if nargin < 11
                  b_L=[];
                  if nargin < 10
                     A=[];
                     if nargin < 9
                        weightY=[];
                        if nargin < 8
                           weightType=[];
                           if nargin < 7
                              t=[];
                              if nargin < 6
                                 x_0=[];
end, end, end, end, end, end, end, end, end, end, end

[m,n] = size(C);

Prob  = ProbDef(1);
Prob.QP = struct('F',[],'c',[],'B',[],'y',[],'Q',[],'R',[],'E',[], ...
    'Ascale',[], 'DualLimit',[],'UseHot',[],'HotFile',[],'HotFreq',[],'HotN',[]);
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
Prob.f_Low = 0;

if length(y) ~= m
   fprintf('Length of y %d, Columns in C are %d\n',length(y),m);
   error('Illegal length of y');
end
if isempty(t)
   Prob.LS.t  = [1:m]';
else
   Prob.LS.t  = full(t(:));
end

if isempty(weightType)
   weightType=0;
end
if weightType==0
   wLS=[];
elseif weightType==1
   if isempty(y) | sum(y)==0
      wLS=[];
   else
      y       = y(:);
      ix      = find(y(:,1)~=0);
      wLS     = zeros(size(y,1),1);
      wLS(ix) = 1./abs(y(ix,1));
   end
elseif weightType==2
   wLS = weightY;
   if length(wLS) ~= m
      size(wLS)
      error('Illegal size of weightY');
   end
end
% Apply weighting directly on C and y
if ~isempty(wLS)
   if size(wLS,1) > 1 & size(wLS,2) > 1
      C   = wLS*C;
      y   = wLS*y;
   else
      wLS = wLS(:);
      C   = diag(wLS)*C;
      y   = wLS.*y;
   end
end

% Store weighting information. Important not to apply it again.
Prob.LS.weightType = 0;
Prob.LS.weightY    = wLS;

% Save C and y after scaling
Prob.LS.C          = C;
Prob.LS.y          = full(y(:));

Prob.LS.SepAlg     = 0;
Prob.LS.damp       = 0;
Prob.LS.L          = [];

if issparse(C)
   Prob.JacPattern = spones(C);
end

Prob.mNonLin  = 0;
Prob.probType = checkType('lls');
Prob.probFile = 0;

% Set Print Level to 0 as default
Prob.PriLevOpt = 0;

global MAX_x MAX_r % Max number of variables/constraints/resids to print

if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

Prob = tomFiles(Prob, 'ls_f', 'ls_g', 'lls_H', [], [], [], 'lls_r', 'lls_J');

% MODIFICATION LOG
%
% 001012  hkh  Written, based on probAssign
% 011204  hkh  Set Prob.f_Low = 0, Prob.LineParam.fLowBnd=0, Prob.LS.SepAlg=0;
% 011204  hkh  Set JacPattern and ConsPattern if sparse matrices
% 021216  hkh  Add fields damp=0 and L empty in Prob.LS
% 030918  ango Default parameters counted wrong, fixed.
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040414  hkh  Use spones to set sparsity patterns in JacPattern, ConsPattern
% 040429  hkh  Safeguard y in LS.y to be a column vector
% 040607  med  Help updates, extra checks removed
% 040728  med  tomFiles used instead
% 041201  hkh  Added check on number of columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 050117  med  MAX_c removed as global, mlint check done
% 050617  hkh  Modified and corrected help
% 060206  med  Added checks for crossover bounds
% 060818  hkh  No need to set Prob.LineParam.fLowBnd=0;, skip
% 061213  med  Moved most input checks to checkAssign
% 110722  hkh Weight y and C if weightType=1 or 2, and store the results
