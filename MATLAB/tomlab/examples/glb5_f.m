% function f = glb4_f(x, Prob)
%
% Shekel 5 test function for global optimization. 
% A and c info are sent using Prob structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Oct 24, 2000.   Last modified Oct 24, 2005.

function f = glb5_f(x, Prob)

% Shekel 5
f=0;
for i = 1:5
    z = x-Prob.user.A(:,i);

    f = f - 1/(z'*z + Prob.user.c(i) );
end

% MODIFICATION LOG
%
% 050117  med  mlint revision