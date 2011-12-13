% rbb_dc - nonlinear constraint gradient matrix 
%          for Rosenbrocks Banana, Problem RB BANANA
%
% function dc = rbbQG_dc(x, Prob)

function dc = rbbQG_dc(x, Prob)

% One row for each constraint, one column for each variable.

dc = [-2*x(1),-1];