% function n = xnargout(S)
%
% Number of function output arguments Matlab 5.2 and later code

% Kenneth Holmstrom, Tomlab Optimization Inc., E-mail: tomlab@tomopt.com.
% Copyright (c) 1997-2007 by Tomlab Optimization Inc., Sweden. $Release: 6.0.0$
% Written Aug 28, 2004.   Last modified Sep 5, 2007.

function n = xnargout(S)

try
   n = abs(nargout(S));
   return
catch
   try
      n = abs(nargout(func2str(S)));
      return
   catch
      err = sprintf('xnargout(''%s'') failed',S);
      error(err);
   end
end

% MODIFICATION LOG
%
% 040828 ango Wrote file
% 041117 ango try-catch construct is default
% 050801 med  Function name changed to nargout
% 070905 med  func2str removed