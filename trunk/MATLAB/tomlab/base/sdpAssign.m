% sdpAssign provides a direct way of setting up linear semidefinite
% programming (SDP) problems with linear matrix inequality
% constraints (LMI) in the TOMLAB (TQ) format.
%
% The information is put into the Tomlab input problem structure Prob.
%
% Prob = sdpAssign( .... )
%
% It is then possible to solve the SDP problem using PENSDP
%
% Result = tomRun('pensdp',Prob);
%
% sdpAssign can also take a predefined PENSDP problem structure as input
% and create a problem structure to be used by TOMLAB.
%
% See the demonstration example in sdpDemo.m in tomlab\examples
%
% -----------------------------------------------------
%
% The LMI minimization problem is defined as
%
%        min   c' * x.  c and x in R^n
%         x
%
%   subject to
%
%      x_L <=   x  <= x_U     x_L,x_U in R^n
%      b_L <= A*x  <= b_U     b_L,b_U in R^ml, A in R^[mlxn]
%
% LMI: Q0^i + sum(i=1:n) (x_j*Q^i_j) <= 0, i = 1, ..., m
%
% For each matrix constraint i, Q0^i, Q^i_j must have the same
% (quadratic) size (or empty)
%
% Equality equations: Set b_L==b_U
% Fixed variables:    Set x_L==x_U
%
% -----------------------------------------------------
%
% Syntax of sdpAssign:
%
% function Prob = sdpAssign(c, SDP, A, b_L, b_U, x_L, x_U, x_0,
%                           Name, fLowBnd);
%
%
% INPUT (One parameter c must always be given)
%
% c             The vector c in c'x in the objective function OR
%               a structure in the PENSDP Structure Format (PSF)
%
% SDP           Structure array describing the LMI constraints.
%               i = 1, ..., m matrix inequalities
%
%   SDP(i).Q{j} are cell-arrays of p x p matrices (or empty)
%   SDP(i).Qidx are column vectors with the indices of the
%               corresponding matrix in SDP(i).Q (or empty)
%
%               All matrices in the same constraints must be of the
%               same dimensions, but two inequalitíes SDP(i) and
%               SDP(j) may have differentyle sized
%               matrices. Constraint i is:
%
%   sum(j=1..k) (x(SDP(i).Qidx(j)) * SDP(i).Q{j})
%
%               where k is the number of non-zero SDP
%               matrices. x(0) is treated as a constant 1 and is
%               not a variable, thus the constant matrix of the
%               SDP constraint is indexed with a zero (0) in
%               the SDP(i).Qidx column vector.
%
%               Any empty Q are interpreted as being all-zero.
%
% A            The linear constraint matrix
% b_L          The lower bounds for the linear constraints
% b_U          The upper bounds for the linear constraints
%
% x_L          Lower bounds on x
% x_U          Upper bounds on x
%
%              b_L, b_U, x_L, x_U must either be empty or of full length
%
% x_0          Starting point x (may be empty)
%
% Name         The name of the problem (string)
%
% fLowBnd      A lower bound on the function value at optimum. Default -1E300
%              Use [] if not known at all. SDP solvers do not use fLowBnd

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.7.0$
% Written July 2, 1999.  Last modified Dec 13, 2006.

function Prob = sdpAssign(c, SDP, A, b_L, b_U, x_L, x_U, x_0, ...
			  Name, fLowBnd)

if nargin < 10
 fLowBnd = [];
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
                 SDP=[];
                 if nargin < 1
                   error('sdpAssign requires at least one parameter c'); 
end,end, end,end,end,end, end,end,end,end

probType = checkType('sdp');

Prob           =ProbDef(1);
Prob.P        = 1;
Prob.probFile = 0;
Prob.probType = probType;
Prob.Name     = deblank(Name);

if isstruct(c)
   Prob.x_0        = c.x0(:);
   Prob.N          = c.vars;
   Prob.QP.c       = c.fobj;
   Prob.PENOPT.pen = c;
   return
end

Prob.QP.c = full(double(c(:)));
Prob.P    = 1;
n         = length(c);
Prob.N    = n;

if ~isempty(fLowBnd) 
   Prob.f_Low=max(Prob.f_Low,fLowBnd); 
end

Prob = checkAssign(Prob, n, x_0, x_L, x_U, b_L, b_U, A);
Prob.mNonLin = 0;
Prob.PENOPT.SDP = SDP;

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
% 020826 hkh  Modified
% 030524 hkh  Use tomFiles to set Prob.FUNCS
% 040102 hkh  Add definition of fields mLin and mNonLin
% 040607 med  Help updates
% 040728 med  tomFiles used instead
% 041201 hkh  Add check on number of columns in A, should be n
% 041222 med  Checking lengths for x_L, x_U and x_0
% 041214 frhe Changed to handle the new LMI format
% 050111 frhe Removed old LMI check.
% 050117 med  mlint revision
% 050127 frhe nProblem and setupFile removed from argument list
% 060206 med  Added checks for crossover bounds
% 060818 hkh  Set Prob.f_Low=max(Prob.f_Low,fLowBnd); if fLowBnd~=[]
% 060822 med  All vectors set to full and double
% 061213 med  Moved most input checks to checkAssign