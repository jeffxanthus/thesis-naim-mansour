% rbb_c - nonlinear constraint vector for Rosenbrocks Banana, Problem RB BANANA
%
% function c = rbbQG_c(x, Prob)

function c = rbbQG_c(x, Prob)

c = -x(1)^2 - x(2);