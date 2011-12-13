% rbb_c - nonlinear constraint vector for Rosenbrocks Banana, Problem RB BANANA
%
% function c = rbb_c(x, Prob)

function cx = rbb_c(x, Prob)

cx = -x(1)^2 - x(2);
