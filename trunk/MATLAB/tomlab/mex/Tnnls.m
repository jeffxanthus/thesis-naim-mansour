%    Tnnls solver, TOMLAB NNLS solver
% 
%    Tnnls solves a linearly constrained least squares
%    problem with both equality and inequality constraints. 
%     
%    Given dense matrices E and A of respective
%    dimensions mE by N and mA by N, and vectors f and b a of
%    respective lengths ME and MA .  This subroutine solves the
%    linearly constrained least squares problem formulated as:
%
%             min   || A x  - b || subject to
%
%                      E x  = f    Equations to be exactly satisfied
%                        x >= 0    Nonnegativitiy constraints
%
%    In case the equality constraints cannot be satisfied, a
%    generalized inverse solution residual vector length is obtained
%    for f-Ex.  This is the minimal length possible for f-Ex.
%
%    Any number of rows in A and E are permitted.
%
% function ...
%   [x, rNorm, mode, Iter, w] = Tnnls ( A, b, E, f, L, ColScale, D, epsRank, epsBlow)
%
% Input: (at least 2 input parameters needed)
%   A         mA x n dense matrix
%   b         mA x 1 dense vector
%   E         mE x n dense matrix
%   f         mE x 1 dense vector
%   L         Variable 1:L are unconstrained. Variable L+1:N non-negative
%   ColScale  If > 0 column scaling of full matrix [A,b;E,f]. Default 0
%   D         If nonempty, n x 1 dense vector with diagonal scaling of columns
%   epsRank   Rank tolerance in least squares part
%   epsBlow   Blow-up parameter, default value SQRT(eps), must be > eps 
%             The reciprocal of this parameter is used in rejecting solution 
%             components as too large when a variable is first brought into the
%             active set. Too large means that the proposed component times the
%             reciprocal of the parameter is not less than the ratio of the 
%             norms of the right-side vector and the data matrix.
%     
% Output:
%   x         Solution x
%   rNorm     || [Ax-b;Ex-f] ||
%   mode      Exit status:      
%             0  Success
%             1  Max. number of iterations (equal to 3*(N-L)) exceeded. 
%                Nearly all problems should complete in fewer than this
%                number of iterations. An approximate solution and its 
%                corresponding residual vector length are in x and rNorm.
%             2  Usage error occurred. Should not be possible.
%   Iter      Number of iterations
%   w         Dual solution vector, Lagrange multipliers
%             w is only computed if a call with 5 output parameters
%

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
%   Written Oct 31, 2000. Last modified Sep 10, 2009.
%# mex

function ...
  [x, rNorm, mode, Iter, w] = Tnnls( A, b, E, f, L, ColScale, D, epsRank, epsBlow)

help Tnnls;
