% qpconQG_c - nonlinear constraint vector for qpcon quick guide
%
% function c = qpconQG_c(x, Prob)

function c = qpconQG_c(x, Prob)

c = sum(x.^2./(1+([1:length(x)]'-1)/3));