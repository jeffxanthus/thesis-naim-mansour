% SSE=funexpm2(x, Prob); sub function
%
% Send xsub and Fe_sub in Prob structure
%
% Test problem by T.Walton

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 21, 1999.     Last modified Nov 26, 2000.

function SSE=funexpm2(x, Prob)

SSE = sum((Prob.user.xSub-x(2)+x(1)*log(1-Prob.user.FeSub)).^2);