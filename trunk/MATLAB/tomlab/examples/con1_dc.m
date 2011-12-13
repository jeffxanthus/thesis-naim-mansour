% function dc = con1_dc(x, Prob)
%
% con1_dc evaluates the constraint gradients dc(x)

function dc = con1_dc(x, Prob)

% The derivative of the two quadratic constraints are linear functions

% One row for each constraint, one column for each variable 

dc = [-x(2)+1 x(2); -x(1)+1 x(1)]';