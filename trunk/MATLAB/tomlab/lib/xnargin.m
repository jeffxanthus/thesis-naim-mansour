% function n = xnargin(S)
%
% Number of function input arguments.

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2008 by Tomlab Optimization Inc., Sweden. $Release: 6.2.0$
% Written Nov 26, 1998.    Last modified Jun 12, 2008.

function n = xnargin(S)

try
   n = abs(nargin(S));
   return;
catch
   try 
      n = abs(nargin(func2str(S)));
      xnargin_cache.(func2str(S)) = n;
      return
   catch
      if(isa(S,'function_handle')), S=func2str(S); end
      err = sprintf('xnargin(''%s'') failed',S);
      fprintf('The user routine %s is most likely incorrectly written. \n', S);
      fprintf('Make sure that it can be executed by itself before proceeding. \n');
      fprintf('For example make the call %s(Prob.x_0, Prob) before the solver call. \n \n', S);
      error(err);
   end
end

% MODIFICATION LOG
%
% 030131 hkh  Add function_handle treatment
% 040728 ango Alternative version useful for mcc users added, commented
% 041117 ango Made the try-catch code default
% 050803 med  Updated error message to provide useful information
% 070905 med  func2str removed
% 080306 pr   Added global variable xnargin_cache, used by tomXtract
% 080612 med  Removed xnargin_cache