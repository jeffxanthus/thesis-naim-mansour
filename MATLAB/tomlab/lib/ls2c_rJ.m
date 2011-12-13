% ls2c_rJ.m
%
% function [r_k, J_k]=ls2c_rJ(x, Prob, varargin)
%
% ls2c_rJ is used for Curve Fitting problems
%
% ls2c_rJ returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem.
% J_k is always a full matrix
%
% ls2c_rJ calls the user function that returns the
% residual, and possibly the Jacobian matrix.
%
% ls2c_rJ is used when implementing the OPTIM TB 2.0 compatibility interface

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written July 28, 1999. Last modified Aug 14, 2006.

function [r_k, J_k]=ls2c_rJ(xNew, Prob, varargin)

rFunc = Prob.FUNCSX.r;

aF=Prob.FUNCSX.funArgIn;

oF=Prob.FUNCSX.funArgOut;

NumDiff=Prob.NumDiff;

x    = Prob.FUNCSX.x;

x(:) = xNew;

if NumDiff >0 | Prob.ADObj > 0 
   % Gradient should not be computed
   if oF < 0
      XDATA = Prob.LS.t;
      r_k=eval(rFunc);
   else
      if aF > 3
         r_k=feval(rFunc, x, Prob.LS.t, Prob, varargin{:} );
      elseif aF > 2
         r_k=feval(rFunc, x, Prob.LS.t, Prob );
      else
         r_k=feval(rFunc, x, Prob.LS.t );
      end
   end
   J_k=[];
elseif NumDiff <0 | Prob.ADObj <0       
   % Hessian should not be computed
   if oF < 0
      XDATA = Prob.LS.t;
      r_k=eval(fFunc);
      J_k = eval(Prob.FUNCSX.J);
   else
      if aF > 3
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t, Prob, varargin{:} );
      elseif aF > 2
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t, Prob );
      else
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t);
      end
   end
else
   if oF < 0
      XDATA = Prob.LS.t;
      f=eval(rFunc);
      J_k = eval(Prob.FUNCSX.J);
      %if oF == -2
      %   d2r = eval(Prob.FUNCSX.d2r);
      %else
      %   d2r = [];
      %end
   else
      if aF > 3
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t, Prob, varargin{:} );
      elseif aF > 2
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t, Prob );
      else
         [r_k,J_k] = feval(rFunc, x, Prob.LS.t );
      end
      %if aF > 1
      %   [r_k,J_k,d2r] = feval(rFunc, x, Prob.LS.t, varargin{:} );
      %else
      %   [r_k,J_k,d2r] = feval(rFunc, x, Prob.LS.t);
      %end
   end
end

r_k=r_k-Prob.LS.y;

% MODIFICATION LOG:
%
% 030129 hkh  Major revision for v4.0. Change principles
% 040522 hkh  3rd arg Prob, then 4th varargin{:}, already split in call
% 060814 med  FUNCSX used for callbacks instead