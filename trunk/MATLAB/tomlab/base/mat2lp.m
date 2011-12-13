% mat2lp writes an LP file. The file is converted from matrices and
% vectors available in MATLAB.
%
%   USAGE:
%
%     mat2lp(Name,Prob);
%
%   INPUT:
%
%    Name          Name of the MPS file with extension.
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
% Written Aug 08, 2005.   Last modified Aug 08, 2005.

function mat2lp(Name, Prob)

if nargin < 2
   error('mat2mps needs two inputs');
end

if issparse(Prob.A)
   Prob.A = [sparse(Prob.QP.c'); Prob.A];
else
   Prob.A = [full(Prob.QP.c'); Prob.A];
end
mpsmex(4, Name, Prob);

% MODIFICATION LOG
%
% 050808 med Written.