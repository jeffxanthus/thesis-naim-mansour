% ls2_rJ.m
%
% function [r_k, J_k]=ls2_rJ(x, Prob, varargin)
%
% ls2_rJ returns both the residual r(x) and the Jacobian J(x) for a
% nonlinear least squares problem.
% J_k is always a full matrix
% 
% ls2_rJ calls the user function that returns the
% residual, and possibly the user function that returns the Jacobian matrix.
%
% ls2_rJ is used when implementing the OPTIM TB 2.x compatibility interface

% Kenneth Holmstrom, Tomlab Optimization Inc, E-mail: tomlab@tomopt.com
% Copyright (c) 1999-2006 by Tomlab Optimization Inc., $Release: 5.5.0$
% Written July 28, 1999. Last modified Aug 14, 2006.

function [r_k, J_k]=ls2_rJ(x, Prob, varargin)

global n_r n_J NARG

nargin;
r=Prob.FUNCSX.r;

n_r = n_r + 1;
if NARG(7) > 1
   r_k = feval(r, x, Prob, varargin{:});
else
   r_k = feval(r, x);
end

if nargout > 1
   J=Prob.FUNCSX.J;
   if ~isempty(J)
      n_J = n_J + 1;
      if NARG(8) > 1
         J_k = full(feval(J, x, Prob, varargin{:}));
      else
         J_k = full(feval(J, x, Prob));
      end
   else
      J_k=[];
   end
end

% MODIFICATION LOG:
%
% 030111 hkh  Revision for v4.0, use NARG
% 060814 med  FUNCS used for callbacks instead