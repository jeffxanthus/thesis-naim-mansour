% function J=L1QG_J(x,Prob)
%
% Compute the Jacobian for the L1 problem Madsen-Tinglett I

function J=L1QG_J(x,Prob)

J = [x(2)-1,x(1);x(2)^2-1,2*x(1)*x(2);x(2)^3-1,3*x(1)*x(2)^2];