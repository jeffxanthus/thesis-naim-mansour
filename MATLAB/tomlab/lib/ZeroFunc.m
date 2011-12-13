% function [f,g] = ZeroFunc(x, Prob)
%
% Callback function called from e.g. Tfzero
%
% If Prob is a structure:
% Prob.FUNCS.f0   The name of the user function
%                (no safety checks are made)
%
% The user function is either defined as
%
%          function f = Func(x);
%          function f = Func(x, Prob);
% or
%          function [f,g] = Func(x);
%          function [f,g] = Func(x, Prob);
%
% If Prob is NOT a structure, then Prob is assumed to be the name of function

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 2000-2006 by Tomlab Optimization Inc., Sweden. $Release: 5.5.0$
% Written Sept 30, 2000.    Last modified Aug 14, 2006.

function [f,g] = ZeroFunc(x, Prob)

if nargout > 1
   Func = Prob.FUNCS.f0;
   if xnargin(Func) >= 2
      [f,g] = feval(Func, x, Prob);
   else
      [f,g] = feval(Func, x);
   end
else
   if isstruct(Prob)
      % Avoid safety tests, will crash if FUNCS.f0 not defined or not function
      Func = Prob.FUNCS.f0;
      if xnargin(Func) >= 2
         f = feval(Func, x, Prob);
      else
         f = feval(Func, x);
      end
   else
      f = feval(Prob, x);
   end
end

% MODIFICATION LOG:
%
% 000930  hkh  Written
% 041201  hkh  Cleanup
% 060814  med  FUNCS used for callbacks instead