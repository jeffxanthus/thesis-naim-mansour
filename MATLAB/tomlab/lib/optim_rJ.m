% optim_rJ.m
%
% function r_k=optim_rJ(xNew, Prob)
%
% optim_rJ is used to implement the OPT TB 2.x interface
%
% The residual r_k is returned
%
% Dependent on if the Jacobian is returned by the user routine or not
% it calls the user routine with one or two output arguments.
%
% If the Jacobian is computed, it is stored in the global variable
% LS_J (and the corresponding x value in LS_xJ)
%
% Test is also made on how many input parameters to use
%
% optim_rJ is called from the TOMLAB gateway function nlp_r.

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1998-2004 by Tomlab Optimization Inc., $Release: 4.7.0$
% Written July 2, 1999.    Last modified Dec 1, 2003.

function r_k=optim_rJ(xNew, Prob)

global LS_xJ LS_J n_J

rFunc=Prob.OPTTB.r;
aF=Prob.OPTTB.funArgIn;
oF=Prob.OPTTB.funArgOut;
NumDiff=Prob.NumDiff;
x    = Prob.OPTTB.x;
x(:) = xNew;

if NumDiff >0 | Prob.ADObj >0 
   % Gradient should not be computed
   if oF < 0
      r_k=eval(rFunc);
   else
      if aF > 1
         r_k=feval(rFunc, x, Prob.varargin{:} );
      else
         r_k=feval(rFunc, x);
      end
   end
   LS_xJ=[];
   LS_J=[];
elseif NumDiff <0 | Prob.ADObj <0       
   % Hessian should not be computed
   if oF < 0
      r_k=eval(fFunc);
      LS_J = eval(Prob.OPTTB.J);
   else
      if aF > 1
         [r_k,LS_J] = feval(rFunc, x, Prob.varargin{:});
      else
         [r_k,LS_J] = feval(rFunc, x);
      end
   end
   n_J = n_J + 1;
   LS_xJ=x(:);
else
   if oF < 0
      f=eval(rFunc);
      LS_J = eval(Prob.OPTTB.J);
   else
      if aF > 1
         [r_k,LS_J] = feval(rFunc, x, Prob.varargin{:});
      else
         [r_k,LS_J] = feval(rFunc, x);
      end
   end
   n_J = n_J + 1;
   LS_xJ=x(:);
end

% MODIFICATION LOG:
%
% 990702  hkh  Written
% 990728  hkh  Global declarations
% 020929  hkh  Must save NumDiff and ADObj before resetting Prob.Psave
% 030126  hkh  Handle matrix x, function handles, inline etc.
% 031101  hkh  Change AutoDiff to new field ADCons