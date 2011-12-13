% rbb_dc - nonlinear constraint gradient matrix 
%          for Rosenbrocks Banana, Problem RB BANANA
%
% function dc = crbb_dc(x, Prob)

function dc = rbb_dc(x, Prob)

% One row for each constraint, one column for each variable.

dc = [-2*x(1),-1];

