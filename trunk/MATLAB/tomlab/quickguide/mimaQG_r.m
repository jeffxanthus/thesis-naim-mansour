% function r=mimaQG_r(x,Prob)
%
% Compute the residual for the minimax problem Madsen-Tinglett II
%
function r=mimaQG_r(x,Prob)

r(1)=-x(1)*(1-x(2)); 
r(2)=-x(1)*(1-x(2)^2);
r(3)=-x(1)*(1-x(2)^3);

% Double the residual, to eliminate the absolute sign
r=[r(:);-r(:)];

r = r - Prob.LS.y;