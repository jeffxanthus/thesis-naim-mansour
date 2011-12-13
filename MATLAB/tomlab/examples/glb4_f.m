% function f = glb4_f(x, Prob)
%
% Shekel 5 test function for global optimization. 
% A and c info are sent using Prob structure

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Oct 24, 2000.   Last modified Oct 24, 2000.

function f = glb4_f(x, Prob)

% Shekel 5

A = Prob.user.A;
c = Prob.user.c;
f = 0;
for i = 1:5
    z = x-A(:,i);
    f = f - 1/(z'*z + c(i) );
end