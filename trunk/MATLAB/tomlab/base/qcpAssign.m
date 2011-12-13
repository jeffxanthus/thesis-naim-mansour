% qcpAssign is a direct way of setting up a Quadratic Complementarity Problem (QCP)
% in the TOMLAB format.
%
% The information is put into the TOMLAB input problem structure Prob.
%
% Prob = qcpAssign(....)
%
% It is then possible to solve the QCP using the TOMLAB KNITRO solver
% knitroTL with tomRun.
%
% See the file tomlab\examples\qcpQG.m for examples.
%
% -----------------------------------------------------
%
% QCP:
%
%        min    0.5*x'*F*x + c' * x.  x in R^n
%
%        s/t   x_L <=   x  <= x_U
%              b_L <= A x  <= b_U
%
% Equality equations: Set b_L==b_U
% Fixed    variables: Set x_L==x_U
%
% MPEC input defines which constraints are complimentary
%
% -----------------------------------------------------
%
%
% Syntax of qcpAssign:
%
% function Prob = qcpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, MPEC, Name, x_min, x_max)
%
% INPUT
%
% F            The square matrix F in the objective function
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
% MPEC         Each row of mpec is one pair. The example below says
%              x(3) _|_ c(2) should be complementary.
%              Exactly two nonzeros per row is allowed.
%
%               mpec = [ ...
%                  3,0,2,0,0,0; ...
%                  4,0,3,0,0,0; ...
%                  5,0,4,0,0,0; ...
%                  ];
%
%              mpec = [   var1,var2 , lin1,lin2 , non1,non2  ; ... ];
%
%              For the linear case only the first 4 columns may be nonzero.
%
% Name         The name of the problem (string)
% x_min        Lower bounds on each x-variable, used for plotting
% x_max        Upper bounds on each x-variable, used for plotting

% Marcus Edvall, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2004-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.7.0$
% Written Apr 18, 2004. Last modified Dec 12, 2006.

function Prob = qcpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, MPEC, Name, x_min, x_max)
                 
if nargin < 12
   x_max=[];
   if nargin < 11
      x_min=[];
      if nargin < 10
         Name=[];
         if nargin < 9
            error('qcpAssign needs at least 9 parameters');
end, end, end, end

global MAX_x MAX_c MAX_r % Max number of variables/constraints/resids to print   

if (~issparse(MPEC))
  MPEC = sparse(MPEC);    % Make sure MPEC is sparse
end

if size(MPEC,2) ~= 6
    error('Input parameter MPEC should have 6 columns');
end

Prob = qpAssign(F, c, A, b_L, b_U, x_L, x_U, x_0, Name,...
                  [], [], [], x_min, x_max);

Prob.probType = checkType('mcp');

Prob = BuildMPEC(Prob,MPEC);

% MODIFICATION LOG
%
% 060524  med  Created
% 061212  med  Help updated (mpec has 6 columns)