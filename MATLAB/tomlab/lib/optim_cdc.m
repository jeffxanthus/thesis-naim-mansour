% optim_cdc.m
%
% function c = optim_cdc(x, Prob)
%
% optim_cdc is used to implement the OPT TB 2.0 interface
%
% The constraint vector c is returned
%
% If the constraint gradient is computed, it is stored in the global variable
% NLP_dc (and the corresponding x value in NLP_xdc)
%
% optim_cdc is called from the TOMLAB gateway function nlp_c.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 29, 1999.  Last modified Dec 22, 2004.

function c = optim_cdc(xNew, Prob)

global NLP_xdc NLP_dc

conFunc = Prob.OPTTB.c;
aF      = Prob.OPTTB.conArgIn;
oF      = Prob.OPTTB.conArgOut;
x       = Prob.OPTTB.x;
x(:)    = xNew;

if Prob.OPTTB.M7Con == 0
   if Prob.ConsDiff > 0 | Prob.ADCons > 0 
      % The Constraint Gradient is to be estimated
      if oF < 0
         ceq  = [];
         c    = eval(conFunc);
      else
         if aF > 1
            [c, ceq]=feval(conFunc, x, Prob.varargin{:});
         else
            [c, ceq]=feval(conFunc, x);
         end
      end

      NLP_xdc=[];
      NLP_dc=[];
   else
      if oF < 0
         ceq  = [];
         dceq = [];
         c    = eval(conFunc);
         dc   = eval(Prob.OPTTB.dc);
      elseif oF >= 4
         if aF > 1
            [c, ceq, dc, dceq] = feval(conFunc, x, Prob.varargin{:});
         else
            [c, ceq, dc, dceq] = feval(conFunc, x);
         end
      else
         dc   = [];
         dceq = [];
         if aF > 1
            [c, ceq] = feval(conFunc, x, Prob.varargin{:});
         else
            [c, ceq] = feval(conFunc, x);
         end
      end
      NLP_xdc=x(:);
      % Keep NLP_dc in the TOMLAB format, one row for each constraint
      if isempty(dceq)
         NLP_dc=dc';
      elseif isempty(dc)
         NLP_dc=dceq';
      else
         NLP_dc=[dceq,dc]';
      end
   end
   c = [ceq(:);c(:)];
else
   if Prob.ConsDiff > 0 | Prob.ADCons > 0 
      % The Constraint Gradient is to be estimated
      if aF > 1
         [c, ceq]=conFunc(x, Prob.varargin{:});
      else
         [c, ceq]=conFunc(x);
      end
      NLP_xdc=[];
      NLP_dc=[];
   else
      if oF >= 4
         if aF > 1
            [c, ceq, dc, dceq] = conFunc(x, Prob.varargin{:});
         else
            [c, ceq, dc, dceq] = conFunc(x);
         end
      else
         dc   = [];
         dceq = [];
         if aF > 1
            [c, ceq] = conFunc(x, Prob.varargin{:});
         else
            [c, ceq] = conFunc(x);
         end
      end
      NLP_xdc=x(:);
      % Keep NLP_dc in the TOMLAB format, one row for each constraint
      if isempty(dceq)
         NLP_dc=dc';
      elseif isempty(dc)
         NLP_dc=dceq';
      else
         NLP_dc=[dceq,dc]';
      end
   end
   c = [ceq(:);c(:)];
end
   
% MODIFICATION LOG:
%
% 030114 hkh Use nargin(conFunc) computed at initialization
% 030115 hkh Add handling of cell, inline and function_handle input
% 030127 hkh Handle matrix x
% 031101 hkh Change AutoDiff to new field ADCons
% 040526 hkh Safe guard c = [ceq(:);c(:)]; if row vectors defined
% 041119 hkh Test oF, allow both two or four outputs from user function
% 041222 hkh Handle new type of function handle in Matlab 7.x