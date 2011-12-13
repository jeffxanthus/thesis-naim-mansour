% minlpQG_dc - nonlinear constraint gradient matrix 
%          for minlp quick guide
%
% function dc = minlpQG_dc(x, Prob)

function dc = minlpQG_dc(x, Prob)

dc = [ ...
         2*x(1)     0.0       1.0  0.0  0.0 ; ...
         0.0   1.5*sqrt(x(2)) 0.0  1.5  0.0];