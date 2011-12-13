% minlpQG_c - nonlinear constraint vector for minlp quick guide
%
% function c = minlpQG_c(x, Prob)

function c = minlpQG_c(x, Prob)

c = [ x(1)^2+x(3) ; sqrt(x(2)^3)+1.5*x(4)];