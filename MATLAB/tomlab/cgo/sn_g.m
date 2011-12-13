% function g = sn_g(y,Prob)
%
% Computes gradient of s_n(y)

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2009 by Tomlab Optimization Inc., $Release: 7.3.0$
% Written Feb 15, 1999.   Last modified June 8, 2008.

function g = sn_g(y,Prob)

if Prob.CGOLIB.usecgolib == 1
  if isnan(Prob.CGO.fnStar)
    g = cgolib(304, Prob.CGOLIB.rbf, y);
  else
    g = 2*(cgolib(304, Prob.CGOLIB.rbf, y) - Prob.CGO.fnStar);    
  end
else
  if isnan(Prob.CGO.fnStar)
    g = tomsol(22,y);
  else
    g = 2*(tomsol(22,y)-Prob.CGO.fnStar);
  end
end

% MODIFICATION LOG
%
% 050427 hkh  Add sign change to f
% 080608 hkh  Remove sign change, slight speedup
% 080616 frhe CGOLIB calls added
