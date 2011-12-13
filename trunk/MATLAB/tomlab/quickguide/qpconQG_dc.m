% qpconQG_dc - nonlinear constraint gradient matrix for qpcon quick guide
%
% function dc = qpconQG_dc(x, Prob)

function dc = qpconQG_dc(x, Prob)

dc = 2*x'./(1+([1:length(x)]-1)/3);