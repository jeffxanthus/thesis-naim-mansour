% function f = glb1_f(x, Prob)
%
% Shekel 5 test function for global optimization. 

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: hkh@tomopt.com.
% Copyright (c) 2000-2004 by Tomlab Optimization Inc., Sweden. $Release: 4.0.0$
% Written Oct 24, 2000.   Last modified Oct 24, 2000.

function f = glb1_f(x, Prob)

% Shekel 5

A = [ 4 4 4 4; 1 1 1 1; 8 8 8 8; 6 6 6 6; 3 7 3 7]';
c = [.1 .2 .2 .4 .4]';
f=0;
for i = 1:5
    z = x-A(:,i);

    f = f - 1/(z'*z + c(i) );
end