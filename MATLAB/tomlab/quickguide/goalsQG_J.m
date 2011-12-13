% function J = goalsQG_J(x, Prob)
%
% Computes the Jacobian to the Multi Criterium Unconstrained & Constrained Nonlinear Problem 
% in point x.
%
% J(i,j) is dr_i/d_x_j

function J = goalsQG_J(x, Prob)

J = [-x(4)  -x(4)  x(4) -x(1)-x(2)+x(3); 
    1-x(2)*x(4)  10+x(4)*x(3)-x(1)*x(4) -1+x(2)*x(4)  1+x(2)*x(3)-x(2)*x(1)];