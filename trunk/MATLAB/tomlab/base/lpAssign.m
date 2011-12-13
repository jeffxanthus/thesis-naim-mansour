% lpAssign is a direct way of setting up a Linear Programming (LP) problem
% in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = lpAssign(...)
%
% It is then possible to solve the LP problem using the TOMLAB LP solver
% lpSimplex with the call:
%    Result = lpSimplex(Prob);
% or
%    Result = tomRun('lpSimplex', Prob, 2);
%
% or any quadratic or general constrained solver:
%
%    Result = tomRun('qpSolve', Prob, 1);
%    Result = tomRun('conSolve', Prob, 1);
%    Result = tomRun('minos', Prob, 1);
%
% It is also possible to run the linprog.m interface, similar to
% linprog in Optimization Toolbox.
%
% See the file tomlab\examples\testlinprog.m for an example
%
% lpAssign may also create an Init File in the TOMLAB Init File format,
% see the input argument setupFile.
% -----------------------------------------------------
%
% LP minimization problem:
%
%
%        min    c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
%
% Syntax of lpAssign:
%
% function Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
%                 setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (One parameter c must always be given)
%
% c            The vector c in c'x in the objective function
% A            The linear constraint matrix
% b_L          The lower bounds for the linear constraints
% b_U          The upper bounds for the linear constraints
% x_L          Lower bounds on x
% x_U          Upper bounds on x
%
%              b_L, b_U, x_L, x_U must either be empty or of full length
%
% x_0          Starting point x (may be empty)
% Name         The name of the problem (string)
% setupFile    The (unique) name as a TOMLAB Init file. If nonempty lpAssign
%              will create a executable m-file with this name and the given
%              problem defined as the first problem in this file.
%              See lp_prob.m, the TOMLAB predefined LP Init File.
%              If empty, no Init File is created.
% nProblem     Number of problems to predefine in the setupFile
%              Not used if setupFile is empty.
%              If empty assumed to be one. Then text are included in the
%              setupFile on how to create additional problems.
% fLowBnd      A lower bound on the function value at optimum. Only used if
%              running nonlinear TOMLAB solvers with line search.
%              Default -1E300
%              A good estimate is not critical. Use [] if not known at all.
% x_min        Lower bounds on each x-variable, used for plotting
% x_max        Upper bounds on each x-variable, used for plotting
% f_opt        Optimal function value(s), if known (Stationary points)
% x_opt        The x-values corresponding to the given f_opt, if known.
%              If only one f_opt, give x_opt as a 1 by n vector
%              If several f_opt values, give x_opt as a length(f_opt) x n matrix
%              If adding one extra column n+1 in x_opt,
%              0 indicates min, 1 saddle (nonlinear problems), 2 indicates max.
%              x_opt and f_opt is used in printouts and plots.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written July 2, 1999. Last modified Jun 6, 2008.

function Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 15
   x_opt=[];
   if nargin < 14
      f_opt=[];
      if nargin < 13
         x_max=[];
         if nargin < 12
            x_min=[];
            if nargin < 11
               fLowBnd=[];
               if nargin < 10
                  nProblem=[];
                  if nargin < 9
                     setupFile=[];
                     if nargin < 8
                        Name=[];
                        if nargin < 7
                           x_0=[];
                           if nargin < 6
                              x_U=[];
                              if nargin < 5
                                 x_L=[];
                                 if nargin < 4
                                    b_U=[];
                                    if nargin < 3
                                       b_L=[];
                                       if nargin < 2
                                          A=[];
                                             if nargin < 1
                           error('lpAssign requires at least one parameter c'); 
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

probType = checkType('lp');
Prob=ProbDef(1);
Prob.QP = struct('F',[],'c',full(double(c(:))),'B',[],'y',[],'Q',[],'R',[],'E',[], ...
    'Ascale',[], 'DualLimit',[],'UseHot',[],'HotFile',[],'HotFreq',[],'HotN',[]);
Prob.P    = 1;
if ~isempty(A)
    if size(A,2) ~= length(c)
        error('Illegal length of c or number of columns in A.');
    end
end
n         = max(size(A,2),length(c));
Prob.N = n;

Prob.probType = probType;
Prob.probFile = 0;

if ~isempty(fLowBnd)
    Prob.f_Low=max(Prob.f_Low,fLowBnd);
end

if length(Prob.QP.c) ~= n
    error('Input c does not have correct length');
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
Prob.mNonLin = 0;

% Set Print Level to 0 as default
Prob.PriLevOpt=0;

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
Prob.Name = deblank(Name);

% MODIFICATION LOG
%
% 990914  hkh  Modify for new automatic setup file generation
% 010528  hkh  Add definition of problem name if no IF file created
% 030524  hkh  Use tomFiles to define the Prob.FUNCS structure
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040607  med  problemName to Name, help fixed for Init File
% 040728  med  tomFiles used instead
% 041123  hkh  Change call to tomRun in help
% 041201  hkh  Add check on columns in A, should be n
% 041222  med  Checking lengths for x_L, x_U and x_0
% 051212  med  More size checks added on inputs
% 060206  med  Added checks for crossover bounds
% 060814  med  FUNCS used for callbacks instead
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822  med  All vectors set to full and double, check for c added
% 061213  med  Moved most input checks to checkAssign
% 080606  med  Removed PrintAssign