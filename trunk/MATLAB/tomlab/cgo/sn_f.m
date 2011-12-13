% function f = sn_f(y,Prob)
%
% Compute s_n(y) using tomsol dll

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Feb 15, 1999.   Last modified June 8, 2008.

function f = sn_f(y, Prob)

if Prob.CGOLIB.usecgolib == 1
  if isnan(Prob.CGO.fnStar)
    f = cgolib(303, Prob.CGOLIB.rbf, y);
  else
    f = (cgolib(303, Prob.CGOLIB.rbf, y)-Prob.CGO.fnStar)^2;
  end
else
  if isnan(Prob.CGO.fnStar)
    f = tomsol(21, y);
  else
    f = (tomsol(21,y)-Prob.CGO.fnStar)^2;
  end
end

% MODIFICATION LOG
%
% 050427 hkh  Add sign change to f
% 080608 hkh  Remove sign change. Slight speedup.
% 080616 frhe CGOLIB calls added.