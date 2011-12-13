% glc5_c.m Compute c(x) for HS 332
%
% function c =glca5_c(x, Prob)

function c = glc5_c(x, Prob)

% HS 332

t = Prob.user.t;
c = max(180/pi*atan(abs((1./t-x(1))./(log(t)+x(2)))));