% SSE=funm2exp(para,xsub,Fe_sub); sub function
%
% Test problem by T.Walton

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1999-2005 by Tomlab Optimization Inc., Sweden. $Release: 4.7.0$
% Written Feb 21, 1999.     Last modified Sep 9, 1999.

function SSE=funm2exp(para,xsub,Fe_sub)

Fi_x0   = -para(1)*log(ones(length(xsub),1)-Fe_sub);
xsubmx0 = xsub-para(2);
SSE     = sum((xsubmx0-Fi_x0).^(2*ones(length(xsub),1)));