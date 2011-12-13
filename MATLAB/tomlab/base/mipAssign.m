% mipAssign is a direct way of setting up a Mixed-Integer Programming (MIP)
% problem in the TOMLAB (TQ) format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = mipAssign(...)
%
% It is then possible to solve the MIP problem using the TOMLAB MIP solver
% mipSolve with the call:  Result = mipSolve(Prob);
%
% or using the cutting plane solver (not as efficient as mipSolve):
%
% Result = cutplane(Prob);
%
% mipAssign may also create an Init File in the TOMLAB Init File format,
% see the input argument setupFile.
%
% -----------------------------------------------------
%
% MIP minimization problem:
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
% x(IntVars) are integer values, IntVars is an index set, a subset of [1:n].
%
% -----------------------------------------------------
%
%
% Syntax of mipAssign:
%
% function Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
%                 setupFile, nProblem, ...
%                 IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
%                 fLowBnd, x_min, x_max, f_opt, x_opt);
%
% INPUT (One parameter c must always be given)
%
% c          The vector c in c'x in the objective function
% A          The linear constraint matrix
% b_L        The lower bounds for the linear constraints
% b_U        The upper bounds for the linear constraints
% x_L        Lower bounds on x
% x_U        Upper bounds on x
%
%            b_L, b_U, x_L, x_U must either be empty or of full length
%
% x_0        Starting point x (may be empty)
% Name       The name of the problem (string)
% setupFile  The (unique) name as a TOMLAB Init File. If nonempty mipAssign
%            will create a executable m-file with this name and the given
%            problem defined as the first problem in this file.
%            See mip_prob.m, the TOMLAB predefined MIP Init File.
%            If empty, no Init File is created.
% nProblem   Number of problems to predefine in the setupFile
%            Not used if setupFile is empty.
%            If empty assumed to be one. Then text are included in the
%            setupFile on how to create additional problems.
%
% --------------------------------------------------------------------------
% The following 5 variables are special for MIP problems, and are assigned to
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
%            A lower value gives higher priority. Setting
%            Prob.MIP.VarWeight = c; for knapsack problems improve convergence.
% KNAPSACK   True if a knapsack problem is to be solved,
%            then a knapsack heuristic is used.
% fIP        An upper bound on the IP value wanted. Makes it possible to
%            cut branches and avoid node computations.
% xIP        The x-values giving the fIP value.
% --------------------------------------------------------------------------
% fLowBnd    A lower bound on the function value at optimum. Default -1E300
%            A good estimate is not critical. Use [] if not known at all.
%            Only used if running some nonlinear TOMLAB solvers with line search
% x_min      Lower bounds on each x-variable, used for plotting
% x_max      Upper bounds on each x-variable, used for plotting
% f_opt      Optimal function value(s), if known (Stationary points)
% x_opt      The x-values corresponding to the given f_opt, if known.
%            If only one f_opt, give x_opt as a 1 by n vector
%            If several f_opt values, give x_opt as a length(f_opt) x n matrix
%            If adding one extra column n+1 in x_opt,
%            0 indicates min, 1 saddle (nonlinear problems), 2 indicates max.
%            x_opt and f_opt is used in printouts and plots.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Nov 15, 1999.   Last modified Oct 23, 2009.

function Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                setupFile, nProblem, ...
                IntVars, VarWeight, KNAPSACK, fIP, xIP, ...
                fLowBnd, x_min, x_max, f_opt, x_opt)

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
                        KNAPSACK=[];
                        if nargin < 12
                           VarWeight=[];
                           if nargin < 11
                              IntVars=[];
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
                         error('mipAssign requires at least one parameter c'); 
end, end, end, end, end, end, end, end, end, end, end, end, end, end, end, ...
end, end, end, end, end 

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print

probType = checkType('mip');

Prob=ProbDef(1);
Prob.QP = struct('F',[],'c',full(double(c(:))),'B',[],'y',[],'Q',[],'R',[],'E',[], ...
    'Ascale',[], 'DualLimit',[],'UseHot',[],'HotFile',[],'HotFreq',[],'HotN',[]);
if ~isempty(A)
   if size(A,2) ~= length(c)
      error('Illegal length of c or number of columns in A.');
   end
end
n      = max(size(A,2),length(c));
Prob.N = n;
Prob   = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);

if any(IntVars == 0) | all(IntVars == 1)
   Prob.MIP = struct('IntVars',logical(abs(IntVars)),...
         'VarWeight',VarWeight,'KNAPSACK',KNAPSACK,'fIP',fIP, 'xIP',xIP,'PI',[],...
         'SC',[],'SI',[],'sos1',[],'sos2',[],'xpcontrol',[],'callback',[]);
else
   Prob.MIP = struct('IntVars',full(double(abs(IntVars(:)))),...
         'VarWeight',VarWeight,'KNAPSACK',KNAPSACK,'fIP',fIP, 'xIP',xIP,'PI',[],...
         'SC',[],'SI',[],'sos1',[],'sos2',[],'xpcontrol',[],'callback',[]);
end

if ~isempty(fLowBnd) 
   Prob.f_Low  = max(Prob.f_Low,fLowBnd); 
end

Prob.probType  = probType;
Prob.probFile  = 0;
Prob.mNonLin   = 0;
Prob.PriLevOpt = 0;
Prob.f_opt     = f_opt;
Prob.x_opt     = x_opt;

if isempty(MAX_x)
   MAX_x=20;
end
if isempty(MAX_c)
   MAX_c=20;
end
if isempty(MAX_r)
   MAX_r=30;
end

Prob           = tomFiles(Prob, 'lp_f', 'lp_g', 'lp_H');
Prob.Name      = deblank(Name);

% MODIFICATION LOG
%
% 000217 hkh  Add definition of special variables
% 010528 hkh  Add definition of problem name if no IF file created
% 011213 hkh  Changed Prob.MIP.SC into SC, SI, semi-continuous and semi-integer
% 011226 hkh  Changed Prob.MIP.SOS1 and SOS2 to sos1 and sos2
% 020822 hkh  Missing apostrophes in save statement
% 021214 hkh  Set probType as lp if isempty(IntVars), otherwise mip
% 030524 hkh  Use tomFiles to define Prob.FUNCS
% 040102 hkh  Add definition of fields mLin and mNonLin
% 040607 med  Help updates, extra checks removed
% 040728 med  tomFiles used instead
% 041201 hkh  Add check on number of columns in A, should be n
% 041222 med  Checking lengths for x_L, x_U and x_0
% 051212 med  More size checks added on inputs
% 060206 med  Added checks for crossover bounds
% 060818 hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822 med  All vectors set to full and double
% 061213 med  Moved most input checks to checkAssign
% 070222 hkh  Revised IntVars input format and definition of Prob.MIP
% 080606 med  Removed PrintAssign
% 091015 hkh  Added definition of Prob.f_opt and Prob.x_opt
% 091023 med  Help updated
