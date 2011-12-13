% qpAssign is a direct way of setting up a Quadratic Programming (QP) problem
% in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = qpAssign(...)
%
% It is then possible to solve the QP problem using the TOMLAB QP solver
% qpSolve with the call:  Result = qpSolve(Prob);
%
% or any general constrained solver:
%
% Result = tomRun('conSolve',Prob);
%
% Adding the parameter 2 gives a call to PrintResult
%
% Result = tomRun('nlpSolve',Prob,2);
%
% If the /SOL toolbox is available, use QPOPT or SQOPT (or SNOPT).
%
% Result = tomRun('qpopt',Prob);
% Result = tomRun('sqopt',Prob);
% Result = tomRun('snopt',Prob);
%
% It is also possible to run the quadprog.m interface, similar to
% quadprog in Optimization Toolbox 2.0.
%
% See the file tomlab\examples\testquadprog.m for an example
%
% qpAssign may also create an Init File in the TOMLAB Init File format,
% see the input argument setupFile.
% -----------------------------------------------------
%
% QP minimization problem:
%
%
%        min   0.5 * x' * F * x + c' * x.  x in R^n
%         x
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of qpAssign:
%
% function Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name,...
%                 setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (One parameter F must always be given. Empty gives default)
%
% F            The matrix F in 0.5 x' F x in the objective function
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
% setupFile    The (unique) name as a TOMLAB Init File. If nonempty qpAssign
%              will create a executable m-file with this name and the given
%              problem defined as the first problem in this file.
%              See qp_prob.m, the TOMLAB predefined QP Init File.
%              If empty, no Init File is created. Also see nProblem.
% nProblem     Number of problems, or problem number, to define in the setupFile
%              Not used if setupFile is empty.
%
%              nProblem = 1 ==> File is created to make it easy to edit new
%              problems into the file. Text are included on how to add new
%              problems. The given problem is set as number 1.
%              If isempty(nProblem) same as nProblem=1.
%
%              length(nProblem) > 1 ==> A file suitable for large test sets
%              are setup, where the problem definition is read from mat-files.
%              Statements for problems nProblem(1) to nProblem(2) are defined.
%              The given input is assumed to be nProblem(1), and the
%              corresponding mat-file is created.
%
%              If nProblem > 1. Additional problems are assumed, and the only
%              thing done is to create a mat-file with the problem.
%
%              If isempty(setupFile), nProblem is not used
%
% fLowBnd      A lower bound on the function value at optimum. Default -1E300
%              A good estimate is not critical. Use [] if not known at all.
%              Only used running some nonlinear TOMLAB solvers with line search
% x_min        Lower bounds on each x-variable, used for plotting
% x_max        Upper bounds on each x-variable, used for plotting
% f_opt        Optimal function value(s), if known (Stationary points)
% x_opt        The x-values corresponding to the given f_opt, if known.
%              If only one f_opt, give x_opt as a 1 by n vector
%              If several f_opt values, give x_opt as a length(f_opt) x n matrix
%              If adding one extra column n+1 in x_opt,
%              0 indicates min, 1 saddle, 2 indicates max.
%              x_opt and f_opt is used in printouts and plots.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written July 2, 1999.   Last modified Jun 6, 2008.

function Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                setupFile, nProblem, fLowBnd, x_min, x_max, f_opt, x_opt)

if nargin < 16
   x_opt=[];
   if nargin < 15
      f_opt=[];
      if nargin < 14
         x_max=[];
         if nargin < 13
            x_min=[];
            if nargin < 12
               fLowBnd=[];
               if nargin < 11
                  nProblem=[];
                  if nargin < 10
                     setupFile=[];
                     if nargin < 9
                        Name=[];
                        if nargin < 8
                           x_0=[];
                           if nargin < 7
                              x_U=[];
                              if nargin < 6
                                 x_L=[];
                                 if nargin < 5
                                    b_U=[];
                                    if nargin < 4
                                       b_L=[];
                                       if nargin < 3
                                          A=[];
                                          if nargin < 2
                                             c=[];
                                             if nargin < 1
                           error('qpAssign requires at least one parameter F'); 
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

probType = checkType('qp');
Prob=ProbDef(1);
Prob.QP = struct('F',F,'c',full(double(c(:))),'B',[],'y',[],'Q',[],'R',[],'E',[], ...
    'Ascale',[], 'DualLimit',[],'UseHot',[],'HotFile',[],'HotFreq',[],'HotN',[]);
if issparse(F)
   Prob.HessPattern = spones(F);
end
Prob.P    = 1;
if ~isempty(c)
   if size(F,1) ~= length(c)
      error('Illegal length of c or number of rows in F.');
   end
end
if ~isempty(A)
   if size(F,1) ~= size(A,2)
      error('Illegal size of F or number of columns in A.');
   end
end

n         = max(size(F,1),length(c));
Prob.N    = n;

if length(Prob.QP.c) ~= n | isempty(Prob.QP.c)
    error('Input c does not have correct length (empty or n)');
end

if size(Prob.QP.F,1) ~= n | size(F,1) ~= size(F,2)
    error('Input F does not have correct size');
end

Prob.probType = probType;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
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

Prob = tomFiles(Prob, 'qp_f', 'qp_g', 'qp_H');
Prob.Name = deblank(Name);

% MODIFICATION LOG
%
% 990908  hkh  Modify for new automatic setup file generation
% 001011  hkh  Set x_0 as zero if empty, otherwise problem in GUI
% 010528  hkh  Add definition of problem name if no IF file created
% 030524  hkh  Use mFiles to define Prob.FUNCS
% 040102  hkh  Add definition of fields mLin and mNonLin
% 040526  hkh  Also set HessPattern if sparse F
% 040607  med  problemName to Name, help fixed for Init File
% 040728  med  tomFiles used instead
% 041123  hkh  Change call to tomRun in help
% 041201  hkh  Add check on number of columns in A, should be n
% 041208  hkh  Move x_min and x_max definition block inside nargin > 0 block
% 041222  med  x_0 safeguard removed, checking lengths for x_L, x_U and x_0
% 050117  med  mlint revision
% 051212  med  More size checks added on inputs
% 060206  med  Added checks for crossover bounds
% 060818  hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822  med  All vectors set to full and double, check for c added
% 061213  med  Moved most input checks to checkAssign
% 080606  med  Removed PrintAssign