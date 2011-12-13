% function f = dace_f(x, Prob)
%
% Likelihood function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Sept 27, 1998. Last modified Sept 13, 2009.

function f = dace_f(x, Prob)

d     = Prob.d;
pEst  = Prob.pEst;
thEst = Prob.thEst;

if pEst == 0
   p = Prob.p;
elseif pEst == 1
   p = x(end)*ones(d,1);
else
   p = x(d+1:d+d);
end

if thEst == 0
   theta = Prob.theta;
elseif thEst == 1
   th    = exp(x(1));
   theta = th*ones(d,1);
else
   theta = exp(x(1:d));
end

ID = Prob.CGOLIB.daceid(Prob.CGOLIB.TRANSFORM+1);

% Call cgolib
cgolib(202, ID, theta, p);

f = cgolib(203, ID);

% MODIFICATION LOG
%
% 040309 hkh  Revise, tomsol(28) returns R(1,1) 1E100 if no inverse possible
% 050218 hkh  Revise, Use SVD to get inv(R) in safe way
% 051101 frhe Old code removed and cgolib call added.
% 060106 frhe Fixed bug with p.
% 060207 ango Change cgolib_mex -> cgolib
% 080112 hkh  Revise for new pEst values, speed up
% 090913 hkh  Revise for new thEst values
