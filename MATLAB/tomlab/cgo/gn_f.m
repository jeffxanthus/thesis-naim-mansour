% function f = gn_f(y, Prob);
%
% Compute gn_f using tomsol dll

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Feb 15, 1999.   Last modified June 8, 2008.

function f = gn_f(y, Prob)

if Prob.CGO.modN == -1
  if(Prob.CGOLIB.usecgolib == 1)
    % Compute mu
    f = cgolib(311, Prob.CGOLIB.rbf, y);
    if f == 0
      f = 1e300;
    else
      f = -1/f;
    end
  else
    % target -inf
    if Prob.CGO.idea == 1
      f = tomsol(26,y,Prob.CGO.fnStar);
    else
      f = tomsol(26,y,Prob.CGO.alpha);
    end
    % Transform with -1/f to avoid the sample points X as stationary points
    if f==0
      f = 1E300;
    else
      f = -1/f;
    end
  end
elseif Prob.CGO.idea == 1
  % fprintf('idea 1\n');
  if Prob.CGOLIB.usecgolib == 1
    f = cgolib(315, Prob.CGOLIB.rbf, y, Prob.CGO.fnStar);
  else
    f = tomsol(23,y, Prob.CGO.fnStar);
  end
else
  % fprintf('idea 2\n');
  if Prob.CGOLIB.usecgolib == 1
    f = cgolib(319, Prob.CGOLIB.rbf, y, Prob.CGO.fnStar);
  else
    f = tomsol(23,y,Prob.CGO.alpha);
  end
end


% MODIFICATION LOG
%
% 050117 med  mlint revision
% 050427 hkh  Add sign change to f
% 050504 hkh  Correction of infStep, both here and in the Fortran code
% 071004 hkh  Transform -inf problem with -1/my, avoid stationary X points
% 080616 frhe CGOLIB calls added