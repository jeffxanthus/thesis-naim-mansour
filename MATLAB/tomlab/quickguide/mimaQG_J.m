% function J=mimaQG_J(x,Prob)
%
% Compute the Jacobian for the minimax problem Madsen-Tinglett II

function J=mimaQG_J(x,Prob)

J = [x(2)-1,x(1);x(2)^2-1,2*x(1)*x(2);x(2)^3-1,3*x(1)*x(2)^2];

% Copy the Jacobian with opposite sign, to eliminate the absolute sign

J = [J;-J];