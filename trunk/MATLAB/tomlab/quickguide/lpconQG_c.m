% lpconQG_c - nonlinear constraint vector for lpcon quick guide
%
% function c = lpconQG_c(x, Prob)

function c = lpconQG_c(x, Prob)

% Circle-Triangle
if x(1)==0 & x(2)==0
    c=1;
else
    c=(x(1)*x(3)+x(2)*x(4))^2/(x(1)^2+x(2)^2) - x(3)^2 - x(4)^2;
end