% minlpQG_H - Hessian matrix for minlp quick guide
%
% function H = minlpQG_H(x, Prob)

function H = minlpQG_H(x, Prob)

H = spalloc(5,5,0);