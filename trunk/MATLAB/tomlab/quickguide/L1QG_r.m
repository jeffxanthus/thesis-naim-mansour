% function r=L1QG_r(x,Prob)
%
% Compute the residual for the L1 problem Madsen-Tinglett I

function r=L1QG_r(x,Prob)

r(1)=-x(1)*(1-x(2)); 
r(2)=-x(1)*(1-x(2)^2);
r(3)=-x(1)*(1-x(2)^3);

r=r(:);

r = r - Prob.LS.y;