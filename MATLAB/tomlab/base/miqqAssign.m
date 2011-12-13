% miqqAssign is a direct way of setting up a
% Mixed-Integer Quadratic Programming (MIQQ) problems with quadratic
% constraints in the TOMLAB format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = miqqAssign( ... )
%
% It is then possible to solve the MIQQ problem using the TOMLAB
% /CPLEX MIQQ solver with the call:
%
% Result = tomRun('cplex',Prob,...);
%
% Adding the two parameters [] and 2 gives a call to PrintResult
%
% Result = tomRun('cplex',Prob),[],2;
%
% miqqAssign may also create an Init File in the TOMLAB Init File format,
% see the input argument setupFile
%
% -----------------------------------------------------
%
% MIQQ minimization problem:
%
%
%        min   0.5 * x' * F * x + c' * x.  x in R^n
%         x
%        s/t   x_L   <=   x    <= x_U
%              b_L   <= A*x    <= b_U
%              x'*Qj*x + aj'*x <= rj   , j=1,2,...,nq
%
%              x(i) in Z for i in I
%
% The aj vector may be zero for any given QC, but the Qj matrix must have at
% least one nonzero element.
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of miqqAssign:
%
% function Prob = miqqAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, qc, ...
%                            IntVars, VarWeight, fIP, xIP, ...
%                            Name, setupFile, nProblem, ...
%                            x_min, x_max, f_opt, x_opt);
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
%
% qc           An arbitrary number of qc's is set.
%
%               qc(1).Q   = sparse( <quadratic coefficient nxn matrix> );
%               qc(1).a   = full  ( <linear coefficient nx1 vector   > );
%               qc(1).r_U = <scalar upper bound> ;
%
%              And similarly for qc(2), ... , qc(n_qc).
%
%              The standard interpretation is x'*Q*x + c'*x <= r_U, but it is
%              possible to define an alternative sense x'*Q*x + c'*x >= r_L
%              by setting qc(i).sense to a nonzero value and specifying a
%              lower bound in qc(i).r_L.
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
% VarWeight    Priorities for each variable in the variable selection phase
%              A lower value gives higher priority.
%
% fIP          An upper bound on the IP value wanted. Makes it possible to
%              cut branches and avoid node computations.
% xIP          The x-values giving the fIP value.
% --------------------------------------------------------------------------
% Name         The name of the problem (string)
%
% setupFile    The (unique) name as a TOMLAB Init File. If nonempty miqpAssign
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
% x_min        Lower bounds on each x-variable, used for plotting
% x_max        Upper bounds on each x-variable, used for plotting
% f_opt        Optimal function value(s), if known (Stationary points)
% x_opt        The x-values corresponding to the given f_opt, if known.
%              If only one f_opt, give x_opt as a 1 by n vector
%              If several f_opt values, give x_opt as a length(f_opt) x n matrix
%              If adding one extra column n+1 in x_opt,
%              0 indicates min, 1 saddle, 2 indicates max.
%              x_opt and f_opt is used in printouts and plots.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written May 21, 2004.   Last modified Oct 23, 2009.

function Prob = miqqAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, qc, ...
                           IntVars, VarWeight, fIP, xIP, ...
                           Name, setupFile, nProblem, ...
                           x_min, x_max, f_opt, x_opt)

if nargin < 20
  x_opt=[];
  if nargin < 19
    f_opt=[];
    if nargin < 18
      x_max=[];
      if nargin < 17
        x_min=[];
        if nargin < 16
          nProblem=[];
          if nargin < 15
            setupFile=[];
            if nargin < 14
              Name=[];
              if nargin < 13
                xIP = [];
                if nargin < 12
                  fIP = [];
                  if nargin < 11
                    VarWeight = [];
                    if nargin < 10
                      IntVars = [];
                      if nargin < 9
                        qc = [];
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
error('miqqAssign requires at least one parameter F'); 
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

probType = checkType('miqq');
Prob=ProbDef(1);
Prob.QP = struct('F',F,'c',full(double(c(:))),'B',[],'y',[],'Q',[],'R',[],'E',[], ...
    'Ascale',[], 'DualLimit',[],'UseHot',[],'HotFile',[],'HotFreq',[],'HotN',[]);
nQC = size(qc,2);
for i=1:nQC
    qc(i).a = qc(i).a(:);
end

Prob.QP.qc = qc;
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
if size(F,1) ~= size(F,2)
   error('Illegal size of F, number of rows should match columns.');
end
n         = max(size(F,1),length(c));
Prob.N    = n;
Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
Prob.probType = probType;
Prob.probFile = 0;

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

nQC = size(Prob.QP.qc,2);
Prob.mNonLin = nQC;

for i=1:nQC
    if isfield(Prob.QP.qc, 'sense')
       if Prob.QP.qc(i).sense ~= 0
          Prob.c_L(i,1) = Prob.QP.qc(i).r_L;
          Prob.c_U(i,1) = inf;
       else
          Prob.c_U(i,1) = Prob.QP.qc(i).r_U; 
          Prob.c_L(i,1) = -inf;
       end
    else
       Prob.c_U(i,1) = Prob.QP.qc(i).r_U; 
       Prob.c_L(i,1) = -inf; 
    end
end

if nQC
    Prob = tomFiles(Prob, 'qp_f', 'qp_g', 'qp_H', 'qp_c', 'qp_dc', 'qp_d2c');
else
    Prob = tomFiles(Prob, 'qp_f', 'qp_g', 'qp_H');
end

if any(IntVars == 0) | all(IntVars == 1)
   Prob.MIP = struct('IntVars',logical(abs(IntVars)),...
         'VarWeight',VarWeight,'KNAPSACK',0,'fIP',fIP, 'xIP',xIP,'PI',[],...
         'SC',[],'SI',[],'sos1',[],'sos2',[],'xpcontrol',[],'callback',[]);
else
   Prob.MIP = struct('IntVars',full(double(abs(IntVars(:)))),...
         'VarWeight',VarWeight,'KNAPSACK',0,'fIP',fIP, 'xIP',xIP,'PI',[],...
         'SC',[],'SI',[],'sos1',[],'sos2',[],'xpcontrol',[],'callback',[]);
end

Prob.Name = deblank(Name);

% MODIFICATION LOG
%
% 040517 med  Created
% 040728 med  tomFiles used instead
% 041201 hkh  Add check on number of columns in A, should be n
% 041222 med  Checking lengths for x_L, x_U and x_0
% 041222 med  Removed extra x_0 check
% 050117 med  mlint revision
% 050901 med  Removed DUNDEE struct, now in ProbCheck, ProbDef
% 051212 med  More size checks added on inputs
% 060206 med  Added checks for crossover bounds
% 060316 med  c_L and c_U set correctly
% 060822 med  All vectors set to full and double
% 061213 med  Moved most input checks to checkAssign
% 070222 hkh  Revised IntVars input format and definition of Prob.MIP
% 080522 med  Change tomFiles call if no quadratic constraints
% 080606 med  Removed PrintAssign
% 091023 med  Help updated