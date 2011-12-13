% mps2mat reads an MPS file. The file is converted into matrices and
% vectors available in MATLAB.
%
%   USAGE:
%
%     Prob = mps2mat(Name);
%
%   INPUT:
%
%    Name          Name of the MPS file with extension.
%
%   OUTPUTS:
%
%    Prob          The following fiels must be set:
%      QP.c        The linear term vector.
%
%      A           The constraint matrix.
%      b_L         The lower bounds of the constraints.
%      b_U         The upper bounds of the constraints.
%
%      x_L         The lower box bounds of x.
%      x_U         The upper box bounds of x.
%
%      IntVars     Logical vector describing what variables that are integer
%                  or binary variables. Empty if the problem is not a mixed
%                  integer problem.

% Marcus Edvall, Tomlab Optimization Inc, E-mail: tomlab@tomlab.biz
% Copyright (c) 2005-2005 by Tomlab Optimization Inc., $Release: 4.9.0$
% Written Aug 25, 2005.   Last modified Aug 25, 2005.

function Prob = mps2mat(Name)

if nargin < 1
   error('mps2mat needs one input');
end

[c,A,b_L,b_U,x_L,x_U,IntVars] = mpsmex(1,Name);

if isempty(IntVars)
   Prob = lpAssign(c, A, b_L, b_U, x_L, x_U, [], Name);
else
   Prob = mipAssign(c, A, b_L, b_U, x_L, x_U, [], Name, [], [], IntVars);
end

% MODIFICATION LOG
%
% 050825 med Written.