% function d2c = con1_d2c(x, lam, Prob)
%
% con1_d2c computes the 2nd part of the Hessian to the Lagrangian function,
%    
%           lam' * d2c(x)
%
% in
%
%   L(x,lam) =   f(x) - lam' * c(x)
% d2L(x,lam) = d2f(x) - lam' * d2c(x)
%
% This problem has two quadratic constraints c(x) = [ c_1(x) ; c_2(x) ]:
%
% c = [-1.5 - x(1)*x(2) + x(1) + x(2); 10 + x(1)*x(2)];
%
% The first derivative matrix consists of four linear expresssions
%
% dc = [-x(2)+1     x(2); 
%       -x(1)+1     x(1)];
%
% where the first column is the derivative of the first constraint
%
% The second derivative is then just constant values
%
% The first lam-value, lam(1) corresponds to the first constraint:
%
%                [ d_2 c_1 / dx_1^2      d_2 c1 / dx_1dx_2  ] 
% lam(1)   *
%                [ d_2 c_1 / dx_2dx_1    d_2 c1 / dx_2^2    ] 
%
% in this case
%
%
%                [  0      -1  ] 
% lam(1)   *
%                [ -1       0  ] 
%
% The second lam-value, lam(2) corresponds to the second constraint:
%
%                [ d_2 c_2 / dx_1^2      d_2 c_2 / dx_1dx_2  ] 
% lam(2)   *
%                [ d_2 c_2 / dx_2dx_1    d_2 c_2 / dx_2^2    ] 
%
% in this case
%
%
%                [  0       1  ] 
% lam(1)   *
%                [  1       0  ] 
%
% The input lam is a column vector. To compute the resulting d2c matrix
% then each element of this matrix, d2c_ij is a scalar product between the
% lam vector and the vector created by taking the ij-element in all
% the second derivative matrices of the constraints
%
% In this case we get
%
% d2c =  [        0  [-1 1]*lam ; 
%        [-1 1]*lam           0]; 
%
% where there are only nonzero scalar products in two places, in both cases
% the scalar product [-1 1]*lam 
%
% There are few solvers that can utilize this type of information,
% but it can speed up the solution process and make it more robust for
% difficult nonlinear constraints.
%
% In TOMLAB it is only nlpSolve and minlpBB in TOMLAB /MINLP that use 
% this kind of information. If not using these solvers, this routine does 
% not need to be defined.

function d2c = con1_d2c(x, lam, Prob)


s = [-1 1]*lam; 
d2c=[ 0  s;  s  0]; 