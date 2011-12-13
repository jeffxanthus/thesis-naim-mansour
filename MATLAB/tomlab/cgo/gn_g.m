% function g = gn_g(y, Prob);
%
% Compute gn_g using cgolib

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2008 by Tomlab Optimization Inc., $Release: 6.2.0$
% Written Feb 15, 1999.   Last modified Nov 5, 2008.

function g = gn_g(y, Prob)

if Prob.CGO.modN == -1
  if(Prob.CGOLIB.usecgolib == 1)
    % Compute mu_g
    % Compute mu (we have a transform of mu and hence need mu as well)
    mu_g = cgolib(312, Prob.CGOLIB.rbf, y);
    mu_f = cgolib(311, Prob.CGOLIB.rbf, y);
    if mu_f == 0
      g = ones(length(y), 1)*1e10;
      %g = ones(length(y), 1)*1e300;
    else
      g = (1/mu_f)^2*mu_g;
    end
  else
    error('gn_g not supported when using tomsol');
  end
elseif Prob.CGO.idea == 1
  % fprintf('idea 1\n');
  if Prob.CGOLIB.usecgolib == 1
    g = cgolib(316, Prob.CGOLIB.rbf, y, Prob.CGO.fnStar);
  else
    error('gn_g not supported when using tomsol');
  end
else
  % fprintf('idea 2\n');
  if Prob.CGOLIB.usecgolib == 1
    error('not implemented yet');
  else
    error('gn_g not supported when using tomsol');
  end
end
if any(isnan(g))
   if Prob.PriLevOpt > 0
      fprintf('gn_g: modN %d Idea %d. ',Prob.CGO.modN,Prob.CGO.idea);
      fprintf('%d NaN elements\n',sum(isnan(g)));
   end
   k = sum(isnan(g));
   g(isnan(g)) = 100*rand(k,1);
   %g(isnan(g)) = 1E10;
   %keyboard
end


% MODIFICATION LOG
%
% 080618 frhe Function written (copied from gn_f.m)
% 081105  hkh  Print NaN info if Prob.PriLevOpt > 0
